"""
This script should take a single cluster_representative as input (e.g. Dusseldorf.2372.BMH_01044)
and the path to the output directory from the HardCORE pipeline.

The script will then find everything it needs by navigating to core_genome/core_genome_DB.cluster_sequences.faa and
core_genome/core_genome_DB.cluster.tsv. It will then create a multifasta file (named after the input cluster_rep) and
extract all relevant .faa sequences from cluster_sequences.faa. Maybe also grab sequences of interest from the .ffn in
the Prokka folder.
"""
import logging
import os
from pathlib import Path

import click

from centreseq.bin.core.accessories import measure
from centreseq.bin.core.summary import cluster_tsv_to_df

script = os.path.basename(__file__)
logger = logging.getLogger()
logging.basicConfig(
    format=f'\033[92m \033[1m {script} %(levelname)s:\033[0m %(message)s ',
    level=logging.INFO)


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    return Path(value)


@click.command(help="Given the path to the HardCORE pipeline root directory and the ID of a cluster representative, "
                    "will create a multi-FASTA containing the sequences for all members of that cluster. Generates "
                    "both an .ffn and .faa file.")
@click.option('-i', '--input_dir',
              type=click.Path(exists=True),
              required=True,
              help='Path to your HardCORE output directory',
              callback=convert_to_path)
@click.option('-o', '--out_dir',
              type=click.Path(exists=False),
              required=True,
              help='Root directory to store all output files',
              callback=convert_to_path)
@click.option('-c', '--cluster_representative',
              type=click.STRING,
              required=False,
              default=None,
              help='Name of the target cluster representative e.g. "Typhi.2299.BMH_00195"')
@measure
def cli(input_dir: Path, out_dir: Path, cluster_representative: str):
    logging.info(f"Extracting member sequences for {cluster_representative}...")

    # Setup
    core_dir = input_dir / 'core_genome'
    core_tsv = core_dir / 'core_genome_DB.cluster.tsv'
    out_faa = out_dir / f"{cluster_representative}_cluster_sequences.faa"
    out_ffn = out_dir / f"{cluster_representative}_cluster_sequences.ffn"

    # Validation
    if validate_folder_contents(input_dir) is False:
        raise Exception("FAILED: Provided HardCORE directory does not match expected structure.")

    # Delete files if they already exist
    if out_faa.exists():
        os.remove(str(out_faa))
    if out_ffn.exists():
        os.remove(str(out_ffn))

    # Make out_dir
    os.makedirs(str(out_dir), exist_ok=True)

    # List of cluster member sequence IDs
    member_list = get_target_member_list(core_tsv, cluster_representative)

    # List of sample names of members
    member_sample_name_list = [member.split("_")[0] for member in member_list]

    # Prepare dict of member_sample_name:member_rep_seq relationship
    member_dict = dict(zip(member_sample_name_list, member_list))

    # Get prokka dir for each
    prokka_dir = input_dir / 'prokka'
    for sample_name, member in member_dict.items():
        sample_prokka_dir = prokka_dir / sample_name
        try:
            faa = list(sample_prokka_dir.glob("*.faa"))[0]
            ffn = list(sample_prokka_dir.glob("*.ffn"))[0]
        except IndexError as e:
            logging.error(f"ERROR: Could not find Prokka files for {sample_name}. Aborting.")
            raise e
        extract_contigs(contigs=faa, target_contig=member, outfile=out_faa)
        extract_contigs(contigs=ffn, target_contig=member, outfile=out_ffn)
    logging.info(f"Done! Results available at {out_dir}")


def validate_folder_contents(root_dir: Path) -> bool:
    prokka_dir = root_dir / 'prokka'
    core_dir = root_dir / 'core_genome'
    mmseqs_dir = root_dir / 'mmseqs2'
    core_tsv = core_dir / 'master_genome_DB.cluster.tsv'
    dir_list = [prokka_dir, core_dir, mmseqs_dir]
    dir_check_list = [d.is_dir() for d in dir_list]
    if False in dir_check_list:
        return False
    if core_tsv.is_file() is False:
        return False
    return True


def get_target_member_list(core_tsv: Path, cluster_representative: str) -> list:
    df = cluster_tsv_to_df(core_tsv)
    df = df[df['cluster_rep'] == cluster_representative]
    member_list = list(df['cluster_member'].values)
    return member_list


def extract_contigs(contigs: Path, target_contig: str, outfile: Path) -> Path:
    outfile_ = open(str(outfile), "a+")
    with open(str(contigs), 'r') as infile:
        write_flag = False
        lines = infile.readlines()
        lines = [line.strip() for line in lines]
        for line in lines:
            if line.startswith(">"):
                if target_contig in line:
                    outfile_.write(line + "\n")
                    write_flag = True
                else:
                    write_flag = False
            elif write_flag:
                outfile_.write(line + "\n")
            else:
                continue
        outfile_.close()
    return outfile


if __name__ == "__main__":
    cli()
