import copy
import csv
import gzip
import logging
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from centreseq.bin.core.accessories import run_subprocess

main_log = logging.getLogger('main_log')


def read_seqs(infile, filter_list=None):
    """
    Reads up sequences from a path to a fasta file
    :param infile: path to fasta file
    :param filter_list: Strings that should be in the description of the sequences
    :return: a list of strings
    """
    r = []
    f = open_possible_gzip(infile)
    for seq in SeqIO.parse(f, "fasta"):
        if filter_list is not None:
            assert isinstance(filter_list, list)
            if any([x in seq.description for x in filter_list]):
                r.append(seq)
        else:
            r.append(seq)
    f.close()
    return r


def faster_fasta_searching(infile, filter_list=[]):
    """
    Loads up a fasta file into a list. Much faster than using SeqIO
    :param infile: fasta infile
    :param filter_list: a list of sequence ids you want to keep. If you want to keep everything pass []
    :return:
    """
    skip = True
    gene_name = ""
    gene_description = ""
    seq = ""
    seqs_all = []
    f = open_possible_gzip(infile)
    for line in f:
        if line[0] == ">":
            # Resolve last gene
            if (filter_list == []) | (gene_name in filter_list):
                seqs_all.append(SeqRecord(Seq(seq), id=gene_name, name=gene_name, description=gene_description))
            # Initialize new gene
            seq = ""
            gene_name = line.split(" ")[0].lstrip(">")
            gene_description = line.rstrip("\n")
            # If we want everything
            if filter_list == []:
                skip = False
            else:
                # Keep this gene
                if gene_name in filter_list:
                    skip = False
                else:
                    skip = True
        elif skip:
            continue
        else:
            # Add sequence to the string
            seq += line.rstrip("\n")
    f.close()

    # Resolve the final gene
    if (filter_list == []) | (gene_name in filter_list):
        seqs_all.append(SeqRecord(Seq(seq), id=gene_name, name=gene_name, description=gene_description))

    return seqs_all


def open_possible_gzip(infile, flags="rt"):
    """
    Opens a file handle for a gzipped or non-zipped file
    :param infile: Path to file
    :param flags:
    :return: file handle
    """
    infile = str(infile)
    if re.search("\.gz$", infile):
        f = gzip.open(infile, flags)
    else:
        f = open(infile, flags)
    return f


def write_seqs_to_file(seq_list, outfile_seq=None):
    """
    Write sequences to file. If not file is given then this is written to a tempfile
    :param seq_list:  a list of sequence objects
    :param outfile_seq: outfile path
    :return: the name of the output file
    """
    if outfile_seq is None:
        outfile_seq = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False).name
    with open(outfile_seq, "w") as f:
        SeqIO.write(seq_list, f, "fasta")
    return outfile_seq


def run_mmseqs(seqs1, seqs2):
    """
    Equivalent to blast_seqs() but uses mmseqs and thus is much faster
    :param seqs1: list of sequences to compare
    :param seqs2: list of sequence to be compared against
    :return:
    """
    query_fasta = write_seqs_to_file(seqs1)
    target_fasta = write_seqs_to_file(seqs2)

    outfile = Path(tempfile.gettempdir()) / (next(tempfile._get_candidate_names()) + ".dat")
    tmpdir = tempfile.TemporaryDirectory()

    # This needs at least mmseqs v8
    result = subprocess.run(["mmseqs"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #m = re.search("MMseqs2 Version: ([0-9])\..+", result.stdout.decode('utf-8'))
    #assert m, "Can't read your mmseqs version, requires at least version 8"
    #assert int(m.group(1)) >= 8, "Require mmseqs at least version 8"

    cmd = f"mmseqs easy-search {query_fasta} {target_fasta} {outfile} {tmpdir.name} --threads 1 --split-memory-limit {max_mem_use} --search-type 3"
    run_subprocess(cmd, get_stdout=True)

    with open(outfile) as f:
        mmseqs_output = f.read().rstrip("\n")

    # I've renamed these for consistency with blast output
    columns = "qseqid,sseqid,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bitscore".split(",")
    return mmseqs_output, columns


def load_pangenome_list(pangenome_list: list):
    """
    Takes a putative core genome list and load it up
    Checks whether there are any paralogs, and therefore, if we need to run the algorithm
    :param pangenome_list: a list of gene names
    """

    # Default is no change
    update = False
    # Check if we have any paralogs, and thus, need to update
    for item in pangenome_list:
        if pd.isnull(item):
            continue
        sp = "_".join(item.split("_")[:-1])
        # If we have any paralogs, then we're going to need to update this
        if len(clusters_global[sp][item]) > 2:
            update = True

    return update


def prehash_clusters(indir: Path):
    """
    Loads up the mmseq self-clustering files (aka list of paralogs)
    Fills self.clusters to be species-> dict from gene to list of genes in its cluster
    """
    clusters = {}
    for sp in os.listdir(str(indir / 'mmseqs2')):
        reg_sub = re.sub('\.[^.]*$', '', sp)
        infile_clusters = f"{str(indir)}/mmseqs2/{sp}/{reg_sub}_DB.cluster.tsv"
        if not Path(infile_clusters).is_file():
            infile_clusters = f"{str(indir)}/mmseqs2/{sp}/{reg_sub}.cluster.tsv"
        clusters[sp] = load_clusters(infile_clusters)
    return clusters


def find_medoid(indir, pangenome_list):
    """
    Finds the medoid of the groups of nucleotide sequences
    Picks a best representative sequence for each species and puts the results into self.cluster
    """
    # If we do need to update this
    seqs_all = []  # ffn sequence objects for everything in the distance matrix
    groups = {}  # sp -> groups of genes

    grouped_seqs = []  # a list of the species annotation for every row/col of the distance matrix
    # Load up the nucleotide sequences such that we can cluster
    for item in pangenome_list:
        if pd.isnull(item):
            continue
        sp = "_".join(item.split("_")[:-1])

        # Load up the nucleotide sequences of those
        infile_ffn = Path(indir / 'prokka' / sp / f"{sp}.ffn")
        seqs = faster_fasta_searching(infile_ffn, filter_list=clusters_global[sp][item])
        groups[sp] = seqs
        seqs_all.extend(seqs)
        grouped_seqs.extend([sp for i in range(len(seqs))])

    # Calculate a distance matrix for the combined group of paralogs
    dist_matrix, all_genes = calc_distance_matrix(seqs_all)
    # The index of the medoid. This index is not in pangenome_list, rather in the full list of sequences (seqs_all)
    medoid_id = define_medoid(dist_matrix, grouped_seqs)
    medoid = seqs_all[medoid_id].id
    main_log.debug(f"Medoid is {str(medoid)}")
    best_seqs = choose_best_seqs(dist_matrix, all_genes, groups, medoid_id)

    # Put best_seqs into the spaces in cluster_post
    cluster_post = copy.copy(pangenome_list)
    j = 0
    for i, e in enumerate(cluster_post):
        if not pd.isnull(e):
            cluster_post[i] = best_seqs[j]
            j += 1
    return medoid, cluster_post


def load_clusters(infile: str, genes_oi: [str] = []):
    """
    Loads up the usearch clusters that have been previously defined
    :param infile: path to cluster file
    :param genes_oi: a list of genes that we want to save
    :return: a dict pointing from each gene, to a list of genes that it is associated with
    """

    cluster = defaultdict(list)
    with open(infile) as f:
        for line in csv.reader(f, delimiter="\t"):
            gene1 = line[0]
            gene2 = line[1]

            if len(genes_oi) > 0:
                if gene1 not in genes_oi:
                    continue
            cluster[gene1].append(gene2)
            if len(genes_oi) > 0:
                if gene2 not in genes_oi:
                    continue
            cluster[gene2].append(gene1)
    # Every cluster needs to have itself in it
    for g in genes_oi:
        assert g in cluster[g]
    return cluster


def calc_distance_matrix(all_seqs: list):
    """
    Uses BLAST to calculate pairwise distances between sequences
    :param all_seqs: A list of sequences to write to file
    :return: a symmatrical matrix of distances and a list of the row/colnames
    """
    # align_output, header = blast_seqs(all_seqs, all_seqs, blast_prog="blastn")
    align_output, header = run_mmseqs(all_seqs, all_seqs)
    all_genes = [x.id for x in all_seqs]
    # Take unique sequences
    dist_matrix = np.ones((len(all_genes), len(all_genes))) * 1e6
    # want bitscore as distance
    for row in align_output.split("\n"):
        # Higher bitscore => smaller distance
        dat = row.split("\t")
        dist_matrix[all_genes.index(dat[header.index('qseqid')]), all_genes.index(
            dat[header.index('sseqid')])] = 1. / float(dat[header.index('bitscore')])
    assert len(all_genes) == dist_matrix.shape[0]
    return dist_matrix, all_genes


def define_medoid(dist_matrix, groups):
    """
    Defines the medoid of a distance matrix
    :param dist_matrix: The numpy array (symmetrical) distance matrix
    :param groups: a list of which sequence matches which species
    :return: the index of the medoid
    """

    # Need to take the minimum over every sequence in a species
    distance_sums = []
    groups = np.array(groups)
    for i in range(dist_matrix.shape[0]):
        distances = []
        for species in set(groups):
            distances.append(np.min(dist_matrix[i, groups == species]))
        distance_sums.append(np.sum(distances))
    return np.argmin(distance_sums)


def choose_best_seqs(dist_matrix, all_genes, groups, medoid):
    """
    Defines the medoid of a distance matrix
    :param dist_matrix: The numpy array (symmetrical) distance matrix
    :param groups: sp -> groups of genes
    :param all_genes: A list of the row/colnames of the distance matrix
    :param medoid: the index of the medoid (from define_medoid)
    :return: the index of the medoid
    """
    assert len(all_genes) == dist_matrix.shape[0]

    distance_to_medoid = np.apply_along_axis(lambda x: np.sum(x - dist_matrix[medoid]), 1, dist_matrix)
    genes_use = []
    for sp in groups:
        genes = [all_genes.index(x.id) for x in groups[sp]]
        # Pick the one with the minimum to the medoid
        gene_oi = all_genes[genes[np.argmin(distance_to_medoid[genes])]]
        genes_use.append(gene_oi)

    # Sort the array, keeping the medoid in front
    assert len(genes_use) == len(groups.keys())
    return genes_use


def determine_memory_avail():
    result = subprocess.run(["cat", "/proc/meminfo"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    m = re.search("MemTotal:\s+([0-9]+)\s+kB", result.stdout.decode('utf-8'))
    assert m, "could not read your memory"
    return float(m.group(1))


def pick_best_nucleotide(summary_report: Path, indir: Path, outdir: Path, n_cpu: int = 1):
    """
    Applies PickBestNucleotide() to each row of summary_report.csv
    :param summary_report: path to summary_report.csv
    :param indir: directory of your HardCORE run
    :param outdir: directory where the summary_report_optimized.csv file will be written
    :param n_cpu: Numbers of cpus. You need at least 10GB of memory per CPU for mmseqs
    :return: the path of the output file
    """
    # Read in the original report then delete
    df_preliminary = pd.read_csv(summary_report, sep="\t")
    summary_report.unlink()

    out_report = outdir / 'summary_report.tsv'

    global pick

    def pick(i):
        pg_list = df_preliminary.iloc[i, 4:].tolist()
        medoid, cluster_post = find_medoid(indir, pg_list)
        return i, medoid, cluster_post

    def update_results(r):
        # input comes from pick()
        i, medoid, post = r
        if (medoid is not None) & (post is not None):
            df_preliminary.iloc[i, 4:] = post
            df_preliminary.iloc[i, 1] = medoid
        pbar.update()

    # since we can't keep this in a class while multiprocessing, we need to save the clusters as a global
    global clusters_global
    clusters_global = prehash_clusters(indir=indir)

    # Max memory for each mmseqs process
    global max_mem_use
    max_mem_use = int(determine_memory_avail() / n_cpu)
    assert max_mem_use >= 10000000, "You need at least 10GB memory per process"

    # First determine those genes that need to be changed at all. This is quick and
    to_change = []
    for i in tqdm.tqdm(range(df_preliminary.shape[0]), desc="Finding clusters to improve"):
        pg_list = df_preliminary.iloc[i, 4:].tolist()
        update = load_pangenome_list(pg_list)
        if update:
            to_change.append(i)

    # Multiprocessing change the core genes

    pool = Pool(n_cpu)
    pbar = tqdm.tqdm(total=len(to_change), desc="Changing cluster reps")
    for i in to_change:
        pool.apply_async(pick, args=(i,), callback=update_results)

    pool.close()
    pool.join()
    pbar.close()

    to_change_names = [df_preliminary["cluster"][i] for i in to_change]

    df_preliminary.to_csv(str(out_report), sep="\t", header=True, index=False)
    return out_report, to_change_names


# if __name__ == "__main__":
#     # pick_best_nucleotide(Path('/mnt/Freya/Analysis/GenomeSimilarity/test_genomes_6/summary_report.csv'),
#     Path("/mnt/Freya/Analysis/GenomeSimilarity/test_genomes_6/"),Path("/mnt/Freya/Analysis/GenomeSimilarity/test_genomes_6/"))
#     # qwe
#     indir_ = "/mnt/Freya/Analysis/GenomeSimilarity/1250-Vibrio_HardCORE/"
#     picker_ = PickBestNucleotide(indir=indir_)
#
#     summary_report_ = indir_ + 'summary_report.csv'
#
#     df_preliminary_ = pd.read_csv(summary_report_, sep="\t")
#     for i_ in [190]:
#         print(i_)
#         pg_list_ = df_preliminary_.iloc[i_, 4:].tolist()
#         picker_.load_pangenome_list(pg_list_)
#         if picker_.update:
#             picker_.find_medoid()
#         print(df_preliminary_.iloc[i_, 0])
#         print("Before {}: After: {}".format(picker_.cluster_pre[1], picker_.cluster_post[1]))
#         print(picker_.cluster_post != picker_.cluster_pre)
#         print(picker_.medoid)
#         print(picker_.cluster_post)
