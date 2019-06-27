import logging
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from collections.__init__ import Counter

from centreseq.bin.core.accessories import run_subprocess

main_log = logging.getLogger('main_log')

"""
TODO: Have some sort of automatic file renaming so Prokka doesn't auto fail
"""


@dataclass
class ProkkaObject:
    """ Dataclass to store metadata on Prokka run """

    # Must be instantiated with these values
    faa: Path
    ffn: Path
    tsv: Path

    # Calculated fields that are determined through __post_init__
    unique_product_list: list = None
    total_product_list: list = None
    product_occurence_dictionary: dict = None
    filtered_product_occurence_dictionary: dict = None

    def __post_init__(self):
        """ Populate calculated parameters. This method is automatically called upon instantiation """
        if self.faa is not None and self.ffn is not None and self.tsv is not None:
            self.populate_unique_gene_list()
            self.populate_total_gene_list()
            self.populate_product_occurence_dictionary()
            self.filter_product_occurence_dictionary()

    def populate_unique_gene_list(self):
        """ Creates a unique gene list from all of the annotated genes found by Prokka """
        df = pd.read_csv(self.tsv, delimiter="\t")
        product_list = list(df['product'].unique())
        product_list = [product for product in product_list if type(product) == str]
        self.unique_product_list = set(product_list)

    def populate_total_gene_list(self):
        """ Grabs a list of all genes in listed by Prokka, including duplicates """
        df = pd.read_csv(self.tsv, delimiter="\t")
        product_list = [product for product in list(df['product']) if type(product) == str]
        self.total_product_list = product_list

    def populate_product_occurence_dictionary(self):
        """ Counts the number of times any product occurs in self.total_gene_list and stores in dict """
        self.product_occurence_dictionary = Counter(self.total_product_list)

    def filter_product_occurence_dictionary(self):
        """ Populates a dictionary with all products that occur more than once in self.product_occurence_dictionary """
        filtered_product_occurence_dictionary = {}
        for key, val in self.product_occurence_dictionary.items():
            if val > 1:
                filtered_product_occurence_dictionary[key] = val
        self.filtered_product_occurence_dictionary = filtered_product_occurence_dictionary


def call_prokka(fasta_path: Path, sample_id: str, outdir: Path, n_cpu: int) -> ProkkaObject:
    """ Makes a system call to Prokka, once complete populates a ProkkaObject with relevant data """
    cmd = f"prokka --centre CORE --compliant --kingdom Bacteria " \
        f"--cpus {n_cpu} --prefix {sample_id} --locustag {sample_id} --outdir {outdir} {fasta_path}"
    run_subprocess(cmd, get_stdout=True)
    # cleanup_prokka(prokka_dir=outdir)  # TODO: Turn this on - will remove extraneous Prokka results
    try:
        faa = list(outdir.glob("*.faa"))[0]
        ffn = list(outdir.glob("*.ffn"))[0]
        tsv = list(outdir.glob("*.tsv"))[0]
        prokka_object = ProkkaObject(faa=faa, ffn=ffn, tsv=tsv)
    except IndexError:
        main_log.error("FAILED: Prokka did not generate the expected files.")
        prokka_object = ProkkaObject(faa=None, ffn=None, tsv=None)
        return prokka_object
    return prokka_object


def prokka_obj_from_results_dir(prokka_dir: Path) -> ProkkaObject:
    try:
        faa = list(prokka_dir.glob("*.faa"))[0]
        ffn = list(prokka_dir.glob("*.ffn"))[0]
        tsv = list(prokka_dir.glob("*.tsv"))[0]
        prokka_object = ProkkaObject(faa=faa, ffn=ffn, tsv=tsv)
    except IndexError:
        main_log.error(f"FAILED: Could not retrieve expected Prokka files from directory {prokka_dir}. "
                       f"This may be a result of input filenames that are too long for Prokka to manage. "
                       f"Additionally, contig names must be <= 37 characters.")
        prokka_object = ProkkaObject(faa=None, ffn=None, tsv=None)
        quit()  # TODO: Test if this is the desired behaviour. Probably best to just fail out of the program here...
        return prokka_object
    return prokka_object


def cleanup_prokka(prokka_dir: Path):
    """ Given a prokka directory, deletes everything except the .ffn, .faa and .tsv files"""
    for f in list(prokka_dir.glob("*")):
        if f.suffix != ".faa" and f.suffix != ".ffn" and f.suffix != ".tsv" and f.is_file():
            f.unlink()


def get_fasta_headers(fasta: Path) -> list:
    """ Pulls headers from any fasta file (e.g. lines starting with >) and returns them as a list """
    fasta_headers = []
    with open(str(fasta)) as f:
        for line in f.readlines():
            if line.startswith(">"):
                line = line[1:].strip()
                fasta_headers.append(line)
    return fasta_headers


def get_product_list_from_prokka_tsv(tsv: Path):
    """ Reads the .tsv file generated by Prokka and returns a unique gene list """
    df = pd.read_csv(tsv, delimiter="\t")
    gene_list = df['product'].unique()
    for gene in gene_list:
        print(gene)


# @measure
def get_product_from_prokka_fasta_header(fasta_header: str) -> str:
    """ Grabs the gene portion of a .ffn or .faa fasta header """
    contig, delim, product = fasta_header.partition(" ")
    return product


# @measure
def get_product_count_dict_from_prokka(ffn_file: Path) -> dict:
    """ Reads Prokka fasta file and returns a dict of counts of each product name extracted from the headers """
    product_list = []
    with open(str(ffn_file)) as f:
        for line in f.readlines():
            if line.startswith(">"):
                line = line.strip()
                product = get_product_from_prokka_fasta_header(line)
                product_list.append(product)
    product_count_dict = dict(Counter(product_list))
    return product_count_dict
