import shutil
from pathlib import Path
from dataclasses import dataclass
from centreseq.bin.core.accessories import run_subprocess, concatenate_faa, sort_faa, log_mmseqs_output
from centreseq.bin.core.sample_handling import SampleObject


@dataclass
class MMseqsObject:
    database: Path
    cluster_database: Path
    cluster_seqs_fasta: Path
    cluster_tsv: Path

    fasta: Path = None
    representative_sequences: Path = None
    representative_sequences_fasta: Path = None
    filtered_cluster_fasta: Path = None


def clean_cluster_fasta(cluster_fasta: Path):
    """
    The output from mmseqs is not a true FASTA and it needs to be processed to remove extraneous lines
    This function expects all FASTA headers to contain the product in addition to the contig identifier. If the product
    is not present in the header line, it will be filtered out.
    """
    filtered_outname = cluster_fasta.with_suffix(".filtered.faa")
    with open(str(cluster_fasta), "r") as f:
        lines = f.readlines()
    with open(str(filtered_outname), "w") as f:
        for line in lines:
            if line.startswith(">"):
                if len(line.split(" ")) > 1:
                    f.write(line)
                else:
                    pass
            else:
                f.write(line)
    return filtered_outname


def call_mmseqs_createdb(fasta: Path, outdir: Path) -> Path:
    """ mmseqs createdb system call """
    out_db = outdir / (fasta.with_suffix("").name + "_DB")
    cmd = f"mmseqs createdb {fasta} {out_db}"
    output = run_subprocess(cmd, get_stdout=True)
    log_mmseqs_output(output)
    return out_db


def call_mmseqs_cluster(database: Path, outdir: Path, n_cpu: int, min_seq_id: float) -> Path:
    """ mmseqs cluster system call """
    out_clusterdb = outdir / (database.with_suffix("").name + "_clusterDB")
    tmp_dir = outdir / "tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)
    cmd = f"mmseqs cluster {database} {out_clusterdb} {tmp_dir} " \
        f"--threads {n_cpu} --min-seq-id {min_seq_id} --add-self-matches"
    output = run_subprocess(cmd, get_stdout=True)
    log_mmseqs_output(output)

    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    return out_clusterdb


def call_mmseqs_linclust(database: Path, outdir: Path, n_cpu: int, min_seq_id: float, alignment_mode: int,
                         coverage_length: float) -> Path:
    """ mmseqs linclust system call """
    out_clusterdb = outdir / (database.with_suffix("").name + "_clusterDB")
    tmp_dir = outdir / "tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)

    """
    linclust option explanation:
    --min-seq-id            only list matches above this sequence identity+
    --alignment-mode 3      0: automatic; 1: only score and end_pos; 2: also start_pos and cov; 3: also seq.id
    --cov-mode 1            coverage of query and target; coverage is defined with -c option
    -c 1                    set c to 100% by default
    --sort-results 1        sort by evalue
    --max-iterations 1000   maximum depth of breadth first search in connected component
    """

    cmd = f"mmseqs linclust {database} {out_clusterdb} {tmp_dir} " \
        f"--threads {n_cpu} --min-seq-id {min_seq_id} --add-self-matches --alignment-mode {alignment_mode} " \
        f"--cov-mode 1 -c {coverage_length} --sort-results 1"

    output = run_subprocess(cmd, get_stdout=True)
    log_mmseqs_output(output)

    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    return out_clusterdb


def call_mmseqs_result2repseq(database: Path, cluster_database: Path, outdir: Path, n_cpu: int) -> Path:
    """ mmseqs result2repseq system call """

    # Only grab representative sequence from clusters
    representative_sequences = outdir / database.with_suffix(".representative_sequences").name
    cmd = f"mmseqs result2repseq {database} {cluster_database} {representative_sequences} --threads {n_cpu}"
    output = run_subprocess(cmd, get_stdout=True)
    log_mmseqs_output(output)
    return representative_sequences


def call_mmseqs_result2flat(database: Path, outdir: Path, cluster_seqs: Path, representative_sequences: bool = False,
                            use_fasta_header: bool = False) -> Path:
    """ mmseqs result2flat system call """

    # Produce FASTA
    if representative_sequences:
        cluster_fasta = outdir / database.with_suffix(".representative_sequences.faa").name
    else:
        cluster_fasta = outdir / database.with_suffix(".cluster_sequences.faa").name

    cmd = f"mmseqs result2flat {database} {database} {cluster_seqs} {cluster_fasta} "
    if use_fasta_header:
        cmd += "--use-fasta-header"
    output = run_subprocess(cmd, get_stdout=True)
    log_mmseqs_output(output)
    return cluster_fasta


def call_mmseqs_result2tsv(database: Path, cluster_database: Path, outdir: Path) -> Path:
    """ mmseqs result2tsv system call """

    # Produce .tsv
    cluster_tsv = outdir / database.with_suffix(".cluster.tsv").name
    cmd = f"mmseqs createtsv {database} {database} {cluster_database} {cluster_tsv} --first-seq-as-repr true"
    output = run_subprocess(cmd, get_stdout=True)
    log_mmseqs_output(output)
    return cluster_tsv


def call_mmseqs_createseqfiledb(database: Path, cluster_database: Path, outdir: Path,
                                min_sequences: int = None, max_sequences: int = None) -> Path:
    """ mmseqs createseqfiledb system call """
    cluster_seq = outdir / (database.with_suffix("").name + "_clusterSEQ")
    cmd = f"mmseqs createseqfiledb {database} {cluster_database} {cluster_seq} "
    if min_sequences is not None:
        cmd += f"--min-sequences {min_sequences} "
    if max_sequences is not None:
        cmd += f"--max-sequences {max_sequences} "
    output = run_subprocess(cmd, get_stdout=True)
    log_mmseqs_output(output)
    return cluster_seq


def self_cluster_pipeline(fasta: Path, outdir: Path, n_cpu: int, min_seq_id: float,
                          coverage_length: float) -> MMseqsObject:
    """ Takes a FASTA file and clusters the sequences within it with mmseqs linclust"""

    # 1. Convert FASTA into mmseqs database format
    database = call_mmseqs_createdb(fasta=fasta, outdir=outdir)

    # 2. Cluster the database file according to min_seq_id. Alignment mode manually set to 3.
    cluster_database = call_mmseqs_linclust(database=database, outdir=outdir, min_seq_id=min_seq_id, n_cpu=n_cpu,
                                            alignment_mode=3, coverage_length=coverage_length)

    # 3. Grab representative sequence from clusters
    representative_sequences = call_mmseqs_result2repseq(database=database, cluster_database=cluster_database,
                                                         outdir=outdir, n_cpu=n_cpu)

    # 4. Create fasta of representative sequences
    representative_sequences_fasta = call_mmseqs_result2flat(database=database,
                                                             cluster_seqs=representative_sequences,
                                                             outdir=outdir, use_fasta_header=True,
                                                             representative_sequences=True)

    # 5. Produce pseudo-FASTA for cluster database
    cluster_seqs = call_mmseqs_createseqfiledb(database=database, cluster_database=cluster_database, outdir=outdir)
    cluster_seqs_fasta = call_mmseqs_result2flat(database=database, cluster_seqs=cluster_seqs, outdir=outdir,
                                                 use_fasta_header=False, representative_sequences=False)

    # 6. Produce .tsv for cluster database containing cluster information
    cluster_tsv = call_mmseqs_result2tsv(database=database, cluster_database=cluster_database, outdir=outdir)

    filtered_cluster_fasta = clean_cluster_fasta(cluster_seqs_fasta)

    representative_sequences_fasta = sort_faa(representative_sequences_fasta)

    mmseqs_object = MMseqsObject(
        database=database,
        cluster_database=cluster_database,
        representative_sequences=representative_sequences,
        representative_sequences_fasta=representative_sequences_fasta,
        cluster_seqs_fasta=cluster_seqs_fasta,
        cluster_tsv=cluster_tsv,
        filtered_cluster_fasta=filtered_cluster_fasta
    )
    return mmseqs_object


def get_core_genome(sample_object_list: [SampleObject], outdir: Path, n_cpu: int, min_seq_id: float,
                    coverage_length: float):
    # Create tmp dir
    tmp_dir = outdir / 'tmp_core_genome'
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)

    faa_list = [sample_object.mmseqs_object.representative_sequences_fasta for sample_object in sample_object_list]
    faa_list = sorted(faa_list)

    core_genome_dir = outdir / 'core_genome'
    if core_genome_dir.exists():
        shutil.rmtree(core_genome_dir)
    core_genome_dir.mkdir(parents=True)

    # Concatenate all .faa files into master_genome, remove it if it already exists
    master_faa_out = core_genome_dir / 'master_genome.faa'
    if master_faa_out.exists():
        master_faa_out.unlink()
    master_faa = concatenate_faa(*faa_list, outname=master_faa_out)
    master_faa = sort_faa(master_faa)

    # MMSEQS CALLS
    # 1. Convert FASTA into mmseqs database format
    database = call_mmseqs_createdb(fasta=master_faa, outdir=core_genome_dir)

    # 2. Cluster the database file according to min_seq_id
    cluster_database = call_mmseqs_linclust(database=database, outdir=core_genome_dir, min_seq_id=min_seq_id,
                                            n_cpu=n_cpu, alignment_mode=3, coverage_length=coverage_length)

    # 3. Produce pseudo-FASTA for cluster database
    cluster_seqs = call_mmseqs_createseqfiledb(database=database, cluster_database=cluster_database,
                                               outdir=core_genome_dir)

    cluster_seqs_fasta = call_mmseqs_result2flat(database=database, cluster_seqs=cluster_seqs, outdir=core_genome_dir,
                                                 use_fasta_header=False, representative_sequences=False)

    # 4. Produce .tsv for cluster database containing cluster information
    cluster_tsv = call_mmseqs_result2tsv(database=database, cluster_database=cluster_database, outdir=core_genome_dir)

    # POST MMSEQS
    # Remove tmp dir
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    # Filter the output FASTA
    filtered_cluster_fasta = clean_cluster_fasta(cluster_seqs_fasta)

    mmseqs_object = MMseqsObject(
        database=database,
        cluster_database=cluster_database,
        cluster_seqs_fasta=cluster_seqs_fasta,
        cluster_tsv=cluster_tsv,
        filtered_cluster_fasta=filtered_cluster_fasta,
        fasta=master_faa
    )

    return mmseqs_object
