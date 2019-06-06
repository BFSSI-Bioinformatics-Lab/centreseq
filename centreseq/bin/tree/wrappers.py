from pathlib import Path

from centreseq.bin.core.accessories import run_subprocess


def call_snp_sites(aligned_multifasta: Path, outdir: Path) -> Path:
    """
    Calls snp-sites on an aligned multiFASTA file and produces a VCF file as output.
    Will only generate an output file if variants are detected.

    https://github.com/sanger-pathogens/snp-sites

    :param aligned_multifasta: Path to multi-FASTA containing alignment of a core gene
    :param outdir: Path to desired output directory
    :return: Path to VCF
    """
    outvcf = outdir / aligned_multifasta.with_suffix(".vcf").name
    cmd = f"snp-sites -v -o {outvcf} {aligned_multifasta}"
    err = run_subprocess(cmd, get_stdout=True)
    return outvcf


def call_muscle(infile: Path) -> Path:
    """
    Produces an aligned version of an input FASTA file (overwrites the original)

    https://www.drive5.com/muscle/
    """
    cmd = f"muscle -in {infile} -out {infile} -maxiters 1"
    run_subprocess(cmd, get_stdout=True)
    return infile
