import logging
import shutil
from multiprocessing import Pool
from pathlib import Path

from tqdm import tqdm

from centreseq.bin.core.accessories import extract_sample_id_from_fasta
from centreseq.bin.core.generate_network_chart import generate_network_chart, generate_network_chart_coding_file
from centreseq.bin.core.mmseqs_wrapper import get_core_genome, self_cluster_pipeline
from centreseq.bin.core.pairwise_comparisons import generate_pairwise_gene_match_report
from centreseq.bin.core.pick_best_nucleotide import pick_best_nucleotide
from centreseq.bin.core.prokka_wrapper import prokka_obj_from_results_dir, call_prokka
from centreseq.bin.core.sample_handling import SampleObject
from centreseq.bin.core.summary import generate_summary_report, generate_core_gene_count_dict, \
    generate_core_gene_count_report, generate_roary_gene_count_dict, generate_roary_gene_count_report
from centreseq.bin.core.validate_clusters import filter_core_cluster_tsv


def core_pipeline(fasta_dir: Path, outdir: Path, n_cpu: int, n_cpu_pickbest: int, min_seq_id: float,
                  coverage_length: float, no_optimize: bool, pairwise: bool):
    """ Generates a core genome and reports """
    main_log = logging.getLogger('main_log')
    main_log.debug(f"minimum sequence identity = {min_seq_id}")
    main_log.debug(f"minimum coverage = {coverage_length}")

    # Prepare [SampleObject]
    sample_object_list = prepare_sample_objects_from_dir(fasta_dir=fasta_dir)
    if len(sample_object_list) < 2:
        main_log.error(f"ERROR: centreseq requires at least 2 samples to proceed. "
                       f"Found {len(sample_object_list)} in {fasta_dir}. Quitting.")
        quit()

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

    # Create report dir
    report_dir = outdir / 'reports'
    report_dir.mkdir(exist_ok=True)

    # SUMMARY REPORT
    main_log.info("Generating summary report")
    summary_report, summary_report_filtered = generate_summary_report(cluster_tsv=core_mmseqs_object.cluster_tsv,
                                                                      core_genome=core_mmseqs_object.fasta,
                                                                      outdir=report_dir)
    main_log.info(f"Summary report available at {summary_report}")

    # NOTE: The user can just do this singleton removal step (and more) very easily in Excel.
    # Not sure if it's worth basically duplicating the summary report.
    main_log.debug(f"Summary report with singletons removed available at {summary_report_filtered}")

    # Call pick_best_nucleotide to improve report
    if not no_optimize:
        main_log.info(f"Optimizing representative sequences by nucleotide")
        summary_report, changed_names = pick_best_nucleotide(summary_report=summary_report, indir=outdir,
                                                             outdir=report_dir,
                                                             n_cpu=n_cpu_pickbest)
        main_log.info(f"Optimized report available at {summary_report}")

        # Log which clusters were altered by pick_best_nucleotide()
        changed_names_log = outdir / 'logs' / 'pick_best_nucleotide_change_log.txt'
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
                                                             outdir=report_dir, coverage_length=coverage_length,
                                                             min_seq_id=min_seq_id)
    main_log.info(f"Core gene count report available {core_gene_count_report}")

    # ROARY-STYLE GENE COUNT REPORT
    roary_gene_count_dict = generate_roary_gene_count_dict(summary_report_tsv=summary_report,
                                                           n_samples=len(sample_object_list))
    roary_gene_count_report = generate_roary_gene_count_report(roary_gene_count_dict=roary_gene_count_dict,
                                                               outdir=report_dir,
                                                               coverage_length=coverage_length,
                                                               min_seq_id=min_seq_id)
    main_log.info(f"Roary-style gene count report available {roary_gene_count_report}")

    # PAIRWISE REPORT
    if pairwise:
        main_log.info("Comparing results from all possible sample pairs")
        pairwise_gene_match_report = generate_pairwise_gene_match_report(summary_report=summary_report,
                                                                         outdir=report_dir)
        main_log.info(f"Pairwise comparison report available at {pairwise_gene_match_report}")

        main_log.info(f"Generating network visualization file")
        sample_id_list = [x.sample_id for x in sample_object_list]
        network_coding = generate_network_chart_coding_file(outdir=outdir, sample_id_list=sample_id_list)
        network_chart = generate_network_chart(pairwise_gene_count_report=pairwise_gene_match_report,
                                               network_coding=network_coding,
                                               roary_report=roary_gene_count_report,
                                               outdir=outdir)
        main_log.info(f"To view the network graph, open a bash terminal, cd to this directory ({outdir}),"
                      f" then run the following command:\n\t\tpython -m http.server\nThis will start a server."
                      f" Now you may open {network_chart} in Chrome, Chromium or Firefox.")
        main_log.info(f"Coding file for network graph available at {network_coding}."
                      f"You can change the values in the group_id column to alter the "
                      f"colouration of the network chart.")
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
                                          min_seq_id=min_seq_id, coverage_length=coverage_length)
    sample.mmseqs_object = mmseqs_object
    return sample, iteration


def prepare_sample_objects_from_dir(fasta_dir: Path):
    """
    Given an input directory containing *.fasta files, creates a SampleObject for each and returns as list
    :param fasta_dir: Directory containing all *.fasta files for core analysis
    :return: Sample dictionary containing sample_id:file_path relationships
    """
    file_list = list(fasta_dir.glob("*.fasta"))
    sample_object_list = []
    for f in file_list:
        sample_object = SampleObject(sample_id=extract_sample_id_from_fasta(f), fasta_path=f)
        sample_object_list.append(sample_object)
    return sample_object_list
