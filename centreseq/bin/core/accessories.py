import itertools
import logging
import multiprocessing
from functools import wraps
from pathlib import Path
from subprocess import PIPE, Popen
from time import time

main_log = logging.getLogger('main_log')
mmseqs_log = logging.getLogger('mmseqs_log')


def get_list_union(*args: list) -> set:
    """ Returns union of x lists as a set"""
    return set().union(*args)


def measure(func):
    """ Decorator to measure the execution time of a function """

    @wraps(func)
    def _time_it(*args, **kwargs):
        start = int(round(time() * 1000))
        try:
            return func(*args, **kwargs)
        finally:
            end_ = int(round(time() * 1000)) - start
            print(f"{func.__name__} execution time: {end_ if end_ > 0 else 0} ms")

    return _time_it


def extract_sample_id_from_fasta(fasta: Path):
    """ Takes path to FASTA file and attempts to assign a Sample ID from the filename """
    sample_id = fasta.with_suffix("").name
    sample_id.replace(" ", "-")
    return sample_id


def get_fasta_headers(fasta: Path) -> list:
    """ Pulls headers any fasta file (e.g. lines starting with >) and returns them as a single list """
    fasta_headers = []
    with open(str(fasta)) as f:
        for line in f.readlines():
            if line.startswith(">"):
                line = line.strip()
                fasta_headers.append(line)
    return fasta_headers


def set_cpu_count(n_cpu: int = None) -> int:
    """
    :param n_cpu: Number of CPUs to set. By default, takes all available - 1.
    :return: Number of threads
    """
    if n_cpu is None:
        n_cpu = multiprocessing.cpu_count() - 1
    return n_cpu


def generate_unordered_pairs(sample_list: list):
    """ Creates a list of of all possible unordered pairs of items in a list """
    return list(itertools.combinations(sample_list, 2))


def extract_contigs(fasta: Path, gene_list: list, outfile: Path) -> Path:
    """
    Searches through a FASTA file and extracts only contigs that are in the provided gene list
    :param fasta: Path to FASTA file
    :param gene_list: List of gene names to extract from contig file
    :param outfile: Path to output file
    :return: Path to output FASTA
    """
    outfile_ = open(str(outfile), "a+")

    write_flag = False
    with open(str(fasta), 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                for gene in gene_list:
                    if gene in line:
                        outfile_.write(line)
                        write_flag = True
                        break
                    else:
                        write_flag = False
            elif write_flag:
                outfile_.write(line)
        outfile_.close()
    return outfile


def run_subprocess(cmd: str, get_stdout=False) -> str:
    """
    :param cmd: System cmd
    :param get_stdout: Set this to true to capture stdout. Will grab either stdout or stderr depending on contents.
    :return: out or err depending on contents, otherwise just returns the cmd provided
    """
    if get_stdout:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        out = out.decode().strip()
        err = err.decode().strip()
        if out != "":
            return out
        elif err != "":
            return err
        else:
            return ""
    else:
        p = Popen(cmd, shell=True)
        p.wait()
        return cmd


def log_mmseqs_output(out: str):
    """
    Drops output from run_subrocess into the logging module
    """
    out_lines = out.split("\n")
    for out_line in out_lines:
        if not out_line.strip() == "":
            mmseqs_log.debug(out_line)


@measure
def filter_gene_count_dict(gene_count_dict: dict) -> dict:
    """
    Filters gene_count_dict to only returns set of key/values with more than 1 hit
    """
    new_gene_count_dict = {gene: count for gene, count in gene_count_dict.items() if count > 1}
    return new_gene_count_dict


def concatenate_faa(*args: [Path], outname: Path):
    """ Takes a list of .faa Path objects and concatenates them into a single file"""
    str_arg_list = [str(f) for f in args]
    outname.touch(exist_ok=True)
    for str_arg in str_arg_list:
        cmd = f"cat {str_arg} >> {outname}"
        run_subprocess(cmd)
    return outname


def sort_faa(faa: Path):
    """
    Sorts input .faa file by length, creates a new sorted version, delete original
    # TODO: Implement this myself instead of using seqkit
    # TODO: Should be very easy to do (sequences to dictionary {header:seq} -> sort keys -> write to file)
    """
    outname = faa.with_suffix(".sorted.faa")
    cmd = f"seqkit sort --by-length --reverse {faa} > {outname}"
    run_subprocess(cmd, get_stdout=True)
    faa.unlink()
    return outname
