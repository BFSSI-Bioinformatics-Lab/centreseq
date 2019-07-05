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

from centreseq.bin.core.summary import cluster_tsv_to_df


def extract_rep_multifasta(indir: Path, outdir: Path, cluster_representative: str):
    logging.info(f"Extracting member sequences for {cluster_representative}...")

    # Setup
    core_dir = indir / 'core_genome'
    core_tsv = core_dir / 'master_genome.cluster.tsv'
    out_faa = outdir / f"{cluster_representative}_cluster_sequences.faa"
    out_ffn = outdir / f"{cluster_representative}_cluster_sequences.ffn"

    # Validation
    if validate_folder_contents(indir) is False:
        raise Exception("FAILED: Provided centreseq directory does not match expected structure.")

    # Delete files if they already exist
    if out_faa.exists():
        os.remove(str(out_faa))
    if out_ffn.exists():
        os.remove(str(out_ffn))

    # Make out_dir
    os.makedirs(str(outdir), exist_ok=True)

    # List of cluster member sequence IDs
    member_list = get_target_member_list(core_tsv, cluster_representative)

    # List of sample names of members
    member_sample_name_list = [member.split("_")[0] for member in member_list]

    # Prepare dict of member_sample_name:member_rep_seq relationship
    member_dict = dict(zip(member_sample_name_list, member_list))

    # Get prokka dir for each
    prokka_dir = indir / 'prokka'
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
    logging.info(f"Done! Results available at {outdir}")


def validate_folder_contents(root_dir: Path) -> bool:
    prokka_dir = root_dir / 'prokka'
    core_dir = root_dir / 'core_genome'
    mmseqs_dir = root_dir / 'mmseqs2'
    core_tsv = core_dir / 'master_genome.cluster.tsv'
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
