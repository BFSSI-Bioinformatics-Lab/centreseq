import itertools
import logging
import warnings
import multiprocessing
from functools import wraps
from pathlib import Path
from subprocess import PIPE, Popen
from time import time

with warnings.catch_warnings():
    warnings.simplefilter('ignore', PendingDeprecationWarning)
    from Bio import SeqIO

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
    handle = open(str(fasta), "r")
    parsed = SeqIO.parse(handle, "fasta")
    fasta_headers = [f">{x.description}" for x in parsed]
    return fasta_headers


def set_cpu_count(n_cpu: int = None) -> int:
    """
    :param n_cpu: Number of CPUs to set. By default, takes all available - 1.
    :return: Number of threads
    """
    if n_cpu is None:
        n_cpu = multiprocessing.cpu_count() - 1
    return n_cpu


def generate_unordered_pairs(sample_list: list) -> list:
    """ Creates a list of of all possible unordered pairs of items in a list """
    return list(itertools.combinations(sample_list, 2))


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


def filter_gene_count_dict(gene_count_dict: dict) -> dict:
    """
    Filters gene_count_dict to only returns set of key/values with more than 1 hit
    """
    new_gene_count_dict = {gene: count for gene, count in gene_count_dict.items() if count > 1}
    return new_gene_count_dict


def concatenate_faa(*args: [Path], outname: Path):
    """ Takes a list of Path objects and concatenates them into a single file"""
    str_arg_list = [str(f) for f in args]
    outname.touch(exist_ok=True)
    for str_arg in str_arg_list:
        cmd = f"cat {str_arg} >> {outname}"
        run_subprocess(cmd)
    return outname


def sort_fasta(fasta: Path, remove_original: bool = False) -> Path:
    original_suffix = fasta.suffix
    handle = open(str(fasta), "r")
    parsed = SeqIO.parse(handle, "fasta")
    sorted_list = [f for f in sorted(parsed, key=lambda x: x.id)]
    sorted_ = fasta.with_suffix(f".sorted{original_suffix}")
    with open(str(sorted_), 'w') as f:
        for s in sorted_list:
            f.write(f">{s.description}\n")
            f.write(f"{str(s.seq)}\n")
    if remove_original:
        fasta.unlink()
    handle.close()
    return sorted_
