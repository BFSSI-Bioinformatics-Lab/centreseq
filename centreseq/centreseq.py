import logging
import os
import shutil
from multiprocessing import Pool
from pathlib import Path

import click
from tqdm import tqdm

from centreseq.bin.core.accessories import prepare_sample_objects_from_dir, set_cpu_count, check_dependencies
from centreseq.bin.core.generate_network_chart import generate_network_chart
from centreseq.bin.core.mmseqs_wrapper import self_cluster_pipeline, get_core_genome
from centreseq.bin.core.pairwise_comparisons import generate_pairwise_gene_match_report
from centreseq.bin.core.pick_best_nucleotide import pick_best_nucleotide
from centreseq.bin.core.prokka_wrapper import call_prokka, prokka_obj_from_results_dir
from centreseq.bin.core.sample_handling import SampleObject
from centreseq.bin.core.summary import generate_summary_report, generate_core_gene_count_dict, \
    generate_core_gene_count_report, \
    generate_roary_gene_count_dict, generate_roary_gene_count_report
from centreseq.bin.core.validate_clusters import filter_core_cluster_tsv

__version__ = "0.1.5"
__authors__ = ["Forest Dussault", "Adrian Verster", "Nicholas Petronella"]
__email__ = "forest.dussault@canada.ca"


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(f"Version\t: {__version__}")
    print(f"Authors\t: {', '.join(__authors__)}")
    print(f"Email\t: {__email__}")
    quit()


def setup_logger(logger_name: str, log_file: Path, title: str, level=logging.INFO):
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter(f'\033[92m \033[1m %(asctime)s - {title} %(levelname)s:\033[0m %(message)s ',
                                  "%Y-%m-%d %H:%M:%S")
    file_formatter = logging.Formatter(f'[%(asctime)s - %(levelname)s]\t %(message)s', "%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(file_formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.setLevel(level)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    return Path(value)


@click.command(help="Given an input directory containing assemblies, establishes a core genome")
@click.option('-f', '--fasta-dir',
              type=click.Path(exists=True),
              required=True,
              help='Path to directory containing *.fasta files for input to the core pipeline',
              callback=convert_to_path)
@click.option('-o', '--outdir',
              type=click.Path(),
              required=True,
              help='Root directory to store all output files. If this directory already exists, the pipeline will '
                   'attempt to skip the Prokka step by reading in existing Prokka output, but will overwrite all other '
                   'existing result files.',
              callback=convert_to_path)
@click.option('-n', '--n-cpu',
              type=click.INT,
              required=False,
              default=None,
              help='Number of CPUs to dedicate to parallelizable steps of the pipeline.'
                   'Will take all available CPUs - 1 by default.')
@click.option('--n-cpu-pickbest',
              type=click.INT,
              required=False,
              default=1,
              help="Number of CPUs for pick_best_nucleotide. You need at least 10GB of RAM per CPU.")
@click.option('-m', '--min-seq-id',
              type=click.FLOAT,
              required=False,
              default=0.95,
              help='Sets the mmseqs cluster parameter "--min-seq-id". Defaults to 0.95.')
@click.option('-c', '--coverage-length',
              type=click.FLOAT,
              required=False,
              default=0.95,
              help='Sets the mmseqs cluster coverage parameter "-c" directly. '
                   'Defaults to 0.95, which is the recommended setting.')
@click.option('--no-optimize',
              is_flag=True,
              default=False,
              help='Set this flag to skip the pick-best-nucleotide step. '
                   'Will improve runtime but decrease quality of results.')
@click.option('--pairwise',
              is_flag=True,
              default=False,
              help='Generate pairwise comparisons of all genomes. '
                   'This output file can be used to view a network chart of the core genome.')
@click.option('-v', '--verbose',
              is_flag=True,
              default=False,
              help='Set this flag to enable more verbose logging.')
@click.option('--version',
              help='Use this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def cli(fasta_dir: Path, outdir: Path, n_cpu: int, n_cpu_pickbest: int, min_seq_id: float, coverage_length: float,
        no_optimize: bool, pairwise: bool, verbose: bool):
    # Create outdir
    os.makedirs(str(outdir), exist_ok=True)

    # Setup logging
    log_dir = outdir / 'logs'
    if log_dir.exists():
        shutil.rmtree(log_dir)
    log_dir.mkdir(parents=True)
    main_log_file = log_dir / 'core_pipeline.log'
    if verbose:
        setup_logger(logger_name='main_log', log_file=main_log_file, title="centreseq", level=logging.DEBUG)
        mmseqs_log_file = log_dir / 'mmseqs2.log'
        setup_logger(logger_name='mmseqs_log', log_file=mmseqs_log_file, title="mmseqs2", level=logging.DEBUG)
    else:
        setup_logger(logger_name='main_log', log_file=main_log_file, title="centreseq", level=logging.INFO)
        mmseqs_log_file = log_dir / 'mmseqs2.log'
        setup_logger(logger_name='mmseqs_log', log_file=mmseqs_log_file, title="mmseqs2", level=logging.INFO)

    main_log = logging.getLogger('main_log')

    # Check dependencies
    dependencies_available = check_dependencies()
    if not dependencies_available:
        quit()

    if min_seq_id > 1 or min_seq_id < 0.02:
        main_log.error("Please enter a --min-seq-id value between 0.02 - 1.00")
        quit()

    # Set CPU count
    if n_cpu is None:
        n_cpu = set_cpu_count(n_cpu)

    # Call the pipeline
    core_pipeline(fasta_dir=fasta_dir, outdir=outdir, n_cpu=n_cpu, n_cpu_pickbest=n_cpu_pickbest, min_seq_id=min_seq_id,
                  coverage_length=coverage_length, no_optimize=no_optimize, pairwise=pairwise)


def core_pipeline(fasta_dir: Path, outdir: Path, n_cpu: int, n_cpu_pickbest: int, min_seq_id: float,
                  coverage_length: float, no_optimize: bool, pairwise: bool):
    main_log = logging.getLogger('main_log')
    main_log.debug(f"minimum sequence identity = {min_seq_id}")
    main_log.debug(f"minimum coverage = {coverage_length}")

    # Prepare [SampleObject]
    sample_object_list = prepare_sample_objects_from_dir(fasta_dir=fasta_dir)

    # Variable to feed to Pool for multiprocessing - this is how many samples will be processed at once
    if len(sample_object_list) < round(n_cpu / 2):
        n_processes = len(sample_object_list)
    else:
        n_processes = round(n_cpu / 4)

    # Check if Prokka results already exist
    # TODO: Check the length of the input FASTA filenames; Prokka fails if they are too long
    prokka_dir = outdir / 'prokka'
    if prokka_dir.exists():
        # Gather existing Prokka results and update corresponding SampleObjects
        main_log.info(f"Detected existing Prokka directory at {prokka_dir}. Attempting to gather Prokka results...")
        prokka_errors = False
        for sample in tqdm(sample_object_list, desc="Collecting Prokka results"):
            expected_prokka_dir = outdir / 'prokka' / sample.sample_id
            if expected_prokka_dir.exists():
                sample.prokka_object = prokka_obj_from_results_dir(prokka_dir=expected_prokka_dir)
            else:
                main_log.error(f"Could not find Prokka results for {sample.sample_id} - directory does not exist")
                prokka_errors = True
        if prokka_errors:
            prokka_errors_message = 'with some errors'
        else:
            prokka_errors_message = 'successfully'
        main_log.info(f'Completed gathering existing Prokka results {prokka_errors_message}')
    else:
        # Async call to prokka_pipeline on every sample
        main_log.info(f"Running Prokka on all provided samples ({n_processes} concurrent processes)")
        pbar = tqdm(total=len(sample_object_list), desc="Running Prokka")
        pool = Pool(processes=n_processes)

        def prokka_callback(result: tuple):
            """ Inner method to update the async progress bar """
            sample_, i_ = result
            results[i_] = sample_
            pbar.update()

        results = ['na'] * len(sample_object_list)
        for i, sample_object in enumerate(sample_object_list):
            pool.apply_async(prokka_pipeline, args=(sample_object, outdir, round(n_cpu / n_processes), i),
                             callback=prokka_callback)

        # Close multiprocessing pool
        pool.close()
        pool.join()
        pbar.close()

        sample_object_list = results
        main_log.info(f"Prokka complete. Successfully annotated {len(sample_object_list)} samples.")

    # Filter out samples where Prokka failed for some reason
    sample_object_list = [sample_object for sample_object in sample_object_list if
                          sample_object.prokka_object.faa is not None]
    sample_object_list = sorted(sample_object_list)

    # Parallel call to mmseqs_pipeline on each sample. Clusters highly similar genes within a sample
    # and produces a representative seq FASTA
    main_log.info("Performing self-clustering on each sample with MMSeqs2")
    pool = Pool(processes=n_processes)
    adjusted_min_seq_id = min_seq_id - 0.01

    def mmseqs_callback(result: tuple):
        """ Inner method to update the async progress bar """
        sample_, i_ = result
        results[i_] = sample_
        pbar.update()

    results = ['na'] * len(sample_object_list)
    pbar = tqdm(total=len(sample_object_list), desc="Running MMseqs2")

    for i, sample_object in enumerate(sample_object_list):
        pool.apply_async(mmseqs_self_cluster_pipeline,
                         args=(
                             sample_object, outdir, round(n_cpu / n_processes), adjusted_min_seq_id, coverage_length,
                             i),
                         callback=mmseqs_callback)

    pool.close()
    pool.join()
    pbar.close()

    # Retrieve result, sort list
    sample_object_list = sorted(results)

    # CORE GENOME
    main_log.info("Generating core genome")
    core_mmseqs_object = get_core_genome(sample_object_list, outdir=outdir, n_cpu=n_cpu, min_seq_id=min_seq_id,
                                         coverage_length=coverage_length)

    # REPAIR CORE GENOME CLUSTER TSV
    main_log.info("Repairing any errors in the core genome cluster tsv file")
    core_mmseqs_object.cluster_tsv = filter_core_cluster_tsv(cluster_tsv=core_mmseqs_object.cluster_tsv,
                                                             outdir=core_mmseqs_object.cluster_tsv.parent)

    # SUMMARY REPORT
    main_log.info("Generating summary report")
    summary_report, summary_report_filtered = generate_summary_report(cluster_tsv=core_mmseqs_object.cluster_tsv,
                                                                      core_genome=core_mmseqs_object.fasta,
                                                                      outdir=outdir)
    main_log.info(f"Summary report available at {summary_report}")
    # NOTE: The user can just do this singleton removal step (and more) very easily in Excel.
    # Not sure if it's worth basically duplicating the summary report.
    main_log.debug(f"Summary report with singletons removed available at {summary_report_filtered}")

    # Call pick_best_nucleotide to improve report
    if not no_optimize:
        main_log.info(f"Optimizing representative sequences by nucleotide")
        summary_report, changed_names = pick_best_nucleotide(summary_report=summary_report, indir=outdir, outdir=outdir,
                                                             n_cpu=n_cpu_pickbest)
        main_log.info(f"Optimized report available at {summary_report}")

        # Log which clusters were altered by pick_best_nucleotide()
        changed_names_log = outdir / 'pick_best_nucleotide_change_log.txt'
        with open(str(changed_names_log), 'w') as f:
            f.write(f"The following entries were changed by pick_best_nucleotide:\n")
            for n in changed_names:
                f.write(f"{n}\n")

    # CORE GENE COUNT REPORT
    core_gene_count_dict = generate_core_gene_count_dict(summary_report_tsv=summary_report,
                                                         n_samples=len(sample_object_list))
    main_log.debug(f"# genes shared among 100% of samples:\t{core_gene_count_dict['n_hard_core_genes']}")
    main_log.debug(f"# genes shared in >=90% of samples:\t{core_gene_count_dict['n_medium_core_genes']}")
    main_log.debug(f"# genes shared in >=50% of samples:\t{core_gene_count_dict['n_soft_core_genes']}")

    core_gene_count_report = generate_core_gene_count_report(core_gene_count_dict=core_gene_count_dict,
                                                             outdir=outdir, coverage_length=coverage_length,
                                                             min_seq_id=min_seq_id)
    main_log.info(f"Core gene count report available {core_gene_count_report}")

    # ROARY-STYLE GENE COUNT REPORT
    roary_gene_count_dict = generate_roary_gene_count_dict(summary_report_tsv=summary_report,
                                                           n_samples=len(sample_object_list))
    roary_gene_count_report = generate_roary_gene_count_report(roary_gene_count_dict=roary_gene_count_dict,
                                                               outdir=outdir,
                                                               coverage_length=coverage_length,
                                                               min_seq_id=min_seq_id)
    main_log.info(f"Roary-style gene count report available {roary_gene_count_report}")

    # PAIRWISE REPORT
    if pairwise:
        main_log.info("Comparing results from all possible sample pairs")
        pairwise_gene_match_report = generate_pairwise_gene_match_report(summary_report=summary_report, outdir=outdir)
        main_log.info(f"Pairwise comparison report available at {pairwise_gene_match_report}")

        main_log.info(f"Generating network visualization file")
        network_chart = generate_network_chart(pairwise_gene_count_report=pairwise_gene_match_report,
                                               roary_report=roary_gene_count_report,
                                               outdir=outdir)
        main_log.info(f"To view the network graph, open a bash terminal, cd to this directory ({outdir}),"
                      f" then run the following command:\n\t\tpython -m http.server\nThis will start a server."
                      f" Now you may open {network_chart} in Chrome, Chromium or Firefox.")
        main_log.info("Done!")


def prokka_pipeline(sample: SampleObject, outdir: Path, n_cpu: int, iteration: int):
    """ Run the call_prokka pipeline, and updates corresponding SampleObject instance """
    prokka_dir = outdir / 'prokka' / sample.sample_id
    prokka_object = call_prokka(sample_id=sample.sample_id, fasta_path=sample.fasta_path, outdir=prokka_dir,
                                n_cpu=n_cpu)
    sample.prokka_object = prokka_object
    return sample, iteration


def mmseqs_self_cluster_pipeline(sample: SampleObject, outdir: Path, n_cpu: int, min_seq_id: float,
                                 coverage_length: float, iteration: int):
    """ Run the call_mmseqs_cluster_pipeline and update corresponding MMSeqsObject instance """
    main_log = logging.getLogger('main_log')
    mmseqs_dir = outdir / 'mmseqs2' / sample.sample_id
    if mmseqs_dir.exists():
        main_log.debug(f"Removing previous mmseqs2 results for {sample.sample_id}")
        shutil.rmtree(mmseqs_dir)
    mmseqs_dir.mkdir(exist_ok=False, parents=True)
    main_log.debug(f"Running mmseqs2 linclust on {sample.sample_id}")
    mmseqs_object = self_cluster_pipeline(fasta=sample.prokka_object.faa, outdir=mmseqs_dir, n_cpu=n_cpu,
                                          min_seq_id=min_seq_id, coverage_length=coverage_length, iterations=3)
    sample.mmseqs_object = mmseqs_object
    return sample, iteration


if __name__ == "__main__":
    cli()
