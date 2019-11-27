import logging
import os
import shutil
from pathlib import Path

import click

from centreseq.bin.core.accessories import set_cpu_count
from centreseq.bin.core.pipelines import core_pipeline
from centreseq.bin.helpers.extract_rep_multifasta import extract_rep_multifasta
from centreseq.bin.helpers.extract_subset import extract_subset
from centreseq.bin.tree.tree_pipeline import tree_pipeline
from centreseq.config import check_dependencies

__version__ = "0.3.5"
__authors__ = ["Forest Dussault", "Adrian Verster", "Nicholas Petronella"]
__email__ = "forest.dussault@canada.ca"

ROOT_DIR = Path(__file__).parent


def print_version(ctx, param, value):
    """ TODO: Activate this once there's something to cite """
    if not value or ctx.resilient_parsing:
        return
    print(f"Version\t: {__version__}")
    print(f"Authors\t: {', '.join(__authors__)}")
    print(f"Email\t: {__email__}")
    print(f"Source\t: https://github.com/BFSSI-Bioinformatics-Lab/centreseq")
    quit()


def setup_logger(logger_name: str, log_file: Path, title: str, level=logging.INFO, stream: bool = True,
                 fancy: bool = True):
    logger = logging.getLogger(logger_name)
    if fancy:
        formatter = logging.Formatter(f'\033[92m \033[1m %(asctime)s - {title} %(levelname)s:\033[0m %(message)s ',
                                      "%Y-%m-%d %H:%M:%S")
        file_formatter = logging.Formatter(f'[%(asctime)s - %(levelname)s]\t %(message)s', "%Y-%m-%d %H:%M:%S")
    else:
        formatter = logging.Formatter()
        file_formatter = logging.Formatter()

    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(file_formatter)

    if stream:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    logger.setLevel(level)
    logger.addHandler(file_handler)


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    return Path(value)


# Base group for centreseq CLI
@click.group()
def centreseq():
    pass


@centreseq.command(short_help="Given an input directory containing assemblies, establishes a core genome",
                   help="Given an input directory containing any number of assemblies (.fasta), centreseq core will "
                        "1) annotate the genomes with Prokka, "
                        "2) perform self-clustering on each genome with MMSeqs2 linclust, "
                        "3) concatenate the self-clustered genomes into a single pan-genome, "
                        "4) cluster the pan-genome with MMSeqs2 linclust, establishing a core genome, "
                        "5) generate helpful reports to interrogate your dataset\n"
                        "Note that if specified output directory already exists, centreseq will search for an existing "
                        "Prokka directory and skip this step if possible.")
@click.option('-f', '--fasta-dir',
              type=click.Path(exists=True),
              required=True,
              help='Path to directory containing *.fasta files for input to the core pipeline',
              callback=convert_to_path)
@click.option('-o', '--outdir',
              type=click.Path(),
              required=True,
              help='Root directory to store all output files. If this directory already exists, the pipeline will '
                   'attempt to skip the Prokka step by reading in the existing Prokka output directory, '
                   'but will overwrite all other existing result files.',
              callback=convert_to_path)
@click.option('-n', '--n-cpu',
              type=click.INT,
              required=False,
              default=None,
              help='Number of CPUs to dedicate to parallelizable steps of the pipeline. '
                   'Will take all available CPUs - 1 by default.')
@click.option('--n-cpu-medoid',
              type=click.INT,
              required=False,
              default=2,
              help="Number of CPUs for the representative medoid picking step (if enabled). "
                   "You will need substantial RAM per CPU.")
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
@click.option('--medoid-repseqs',
              is_flag=True,
              default=False,
              help='This setting will identify the representative medoid nucleotide sequence for each core cluster. '
                   'Enabling this will increase computation time considerably. Note that this parameter has no effect'
                   ' on the number of core clusters detected.')
@click.option('--pairwise',
              is_flag=True,
              default=False,
              help='Generate pairwise comparisons of all core genomes. This setting allows for viewing an  '
                   'interactive network chart which visualizes core genome relatedness.')
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
def core(fasta_dir: Path, outdir: Path, n_cpu: int, n_cpu_medoid: int, min_seq_id: float, coverage_length: float,
         medoid_repseqs: bool, pairwise: bool, verbose: bool):
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
        setup_logger(logger_name='mmseqs_log', log_file=(log_dir / 'mmseqs2.log'), title="mmseqs2", level=logging.DEBUG,
                     fancy=False)
        setup_logger(logger_name='prokka_log', log_file=(log_dir / 'prokka.log'), title="prokka", level=logging.DEBUG,
                     fancy=False)
    else:
        setup_logger(logger_name='main_log', log_file=main_log_file, title="centreseq", level=logging.INFO)
        setup_logger(logger_name='mmseqs_log', log_file=(log_dir / 'mmseqs2.log'), title="mmseqs2", level=logging.INFO,
                     stream=False, fancy=False)
        setup_logger(logger_name='prokka_log', log_file=(log_dir / 'prokka.log'), title="prokka", level=logging.INFO,
                     stream=False, fancy=False)

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
    core_pipeline(fasta_dir=fasta_dir, outdir=outdir, n_cpu=n_cpu, n_cpu_pickbest=n_cpu_medoid, min_seq_id=min_seq_id,
                  coverage_length=coverage_length, medoid_repseqs=medoid_repseqs, pairwise=pairwise)


@centreseq.command(short_help="Produces output for phylogenetic tree software",
                   help="Processes centreseq core output files to produce files that can be fed into phylogenetic tree "
                        "building software.")
@click.option('-s', '--summary-report',
              type=click.Path(exists=True),
              required=True,
              help='Path to summary_report.csv file produced by the core pipeline',
              callback=convert_to_path)
@click.option('-p', '--prokka-dir',
              type=click.Path(exists=True),
              required=True,
              help='Path to the Prokka output directory generated by the core pipeline',
              callback=convert_to_path)
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              help='Root directory to store all output files',
              callback=convert_to_path)
@click.option('-pct', '--percentile',
              type=click.FLOAT,
              required=False,
              default=99.0,
              help='Filter summary report by n_members to the top nth percentile. Defaults to 99.0.')
@click.option('-n', '--n-cpu',
              type=click.INT,
              required=False,
              default=None,
              help='Number of CPUs to dedicate to parallelizable steps of the pipeline.'
                   'Will take all available CPUs - 1 if not specified.')
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
def tree(summary_report: Path, prokka_dir: Path, outdir: Path, percentile: float, n_cpu: int, verbose: bool):
    # Outdir validation
    if outdir.exists():
        logging.error(f"ERROR: Directory {outdir} already exists!")
        quit()
    outdir.mkdir()

    # Setup logging
    log_dir = outdir / 'logs'
    log_dir.mkdir(parents=True)
    main_log_file = log_dir / 'tree_pipeline.log'
    if verbose:
        setup_logger(logger_name='main_log', log_file=main_log_file, title="TreePipeline", level=logging.DEBUG)
    else:
        setup_logger(logger_name='main_log', log_file=main_log_file, title="TreePipeline", level=logging.INFO)
    main_log = logging.getLogger('main_log')

    # Set CPU count
    if n_cpu is None:
        n_cpu = set_cpu_count(n_cpu)

    main_log.info("Started Core Tree Pipeline")
    main_log.debug(f"summary_report:\t{summary_report}")
    main_log.debug(f"prokka_dir:\t{prokka_dir}")
    main_log.debug(f"outdir:\t{outdir}")
    main_log.debug(f"percentile:\t{percentile}")
    main_log.debug(f"n_cpu:\t{n_cpu}")
    main_log.debug(f"verbose:\t{verbose}")

    # Call the pipeline
    tree_pipeline(summary_report=summary_report,
                  prokka_dir=prokka_dir,
                  outdir=outdir,
                  n_cpu=n_cpu,
                  percentile=percentile)


@centreseq.command(short_help="Helper tool to extract sequences from a particular core cluster",
                   help="Given the path to the centreseq core directory and the ID of a cluster "
                        "representative, will create a multi-FASTA containing the sequences for all members of that "
                        "cluster. Generates both an .ffn and .faa file.")
@click.option('-i', '--indir',
              type=click.Path(exists=True),
              required=True,
              help='Path to your centreseq output directory',
              callback=convert_to_path)
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              help='Root directory to store all output files',
              callback=convert_to_path)
@click.option('-c', '--cluster_representative',
              type=click.STRING,
              required=True,
              default=None,
              help='Name of the target cluster representative e.g. "Typhi.2299.BMH_00195"')
@click.option('--version',
              help='Use this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def extract(indir, outdir, cluster_representative):
    extract_rep_multifasta(indir=indir, outdir=outdir, cluster_representative=cluster_representative)


@centreseq.command(short_help="Subset summary_report.tsv to only samples of interest",
                   help="Given an input text file of Sample IDs and a summary report, will return a filtered version "
                        "of the summary report for clusters that belong exclusively in the input sample ID list")
@click.option('-i', '--input-samples',
              type=click.Path(exists=True),
              required=True,
              help='Path to a new line separated text file containing each Sample ID to target',
              callback=convert_to_path)
@click.option('-s', '--summary-report',
              type=click.Path(exists=True),
              required=True,
              help='Path to summary report generated by the centreseq core command, i.e. summary_report.tsv',
              callback=convert_to_path)
@click.option('-o', '--outpath',
              type=click.Path(exists=False),
              required=False,
              default=None,
              help='Path to desired output file. If no value is provided, will create a new report in the same '
                   'directory as the input summary report.',
              callback=convert_to_path)
def subset(input_samples: Path, summary_report: Path, outpath: Path):
    extract_subset(input_samples=input_samples, summary_report=summary_report, outpath=outpath)


@click.command(cls=click.CommandCollection,
               sources=[centreseq],
               help="centreseq builds an annotated core genome using assemblies as input.")
@click.option('--version',
              help='Print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              default=False,
              expose_value=False)
# @click.option('--citation',
#               help='Print the citation for this software and exit.',
#               is_flag=True,
#               is_eager=True,
#               callback=print_citation,
#               default=False,
#               expose_value=False)
def cli():
    pass


if __name__ == "__main__":
    cli()
