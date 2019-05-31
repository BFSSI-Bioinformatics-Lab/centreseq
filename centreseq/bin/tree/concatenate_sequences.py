import multiprocessing
from pathlib import Path
from copy import deepcopy


def concatenate_sequence_directory(sample_ids: [str], sequence_directory: Path, n_processes: int, outdir: Path) -> Path:
    """
    Given a sequence directory containing aligned multi-FASTA files, will attempt to concatenate all of the sequences
    into a single file. Minimal example showing input and output (concatenated_sequences.fasta):

    files in directory:
        gene_1.fasta
        gene_2.fasta

    contents of gene1.fasta:
        >sample_1
        ATCG
        >sample2
        ATTT

    contents of gene2.fasta:
        >sample_1
        GGGGAGGGGGTCA

    concatenated_sequences.fasta:
        >sample_1
        ATCGGGGGAGGGGGTCA
        >sample_2
        ATTTNNNNNNNNNNNNN

    - Provided sample_ids list must contain exact matches to headers in FASTA files
    - All sequences in a single FASTA file must be the exact same length; this is done through an aligner like MUSCLE
    - Expectation is that every FASTA file will have the same headers, or a subset of headers present in sample_ids
    """
    concat_seqs_dict = generate_concat_seqs_dict(sample_ids=sample_ids, indir=sequence_directory,
                                                 n_processes=n_processes)
    outfile = write_concat_seqs_dict(concat_seqs_dict=concat_seqs_dict, outdir=outdir)
    return outfile


def write_concat_seqs_dict(concat_seqs_dict: dict, outdir: Path) -> Path:
    outfile = outdir / 'concatenated_sequences.fasta'
    if outfile.exists():
        outfile.unlink()
    outfile_ = open(str(outfile), 'a+')
    for sample_id, sequence in concat_seqs_dict.items():
        outfile_.write(f">{sample_id}\n")
        outfile_.write(f"{sequence}\n")
    outfile_.close()
    return outfile


def generate_concat_seqs_dict(sample_ids: set, indir: Path, n_processes=4) -> dict:
    # Potentially takes a lot of RAM. Stores the concatenated sequences for each sample.
    sequence_storage = {sample_id: "" for sample_id in sample_ids}
    fasta_files = sorted(list(indir.glob("*.fasta")) + list(indir.glob("*.ffn")))

    # Set # of concurrent processes to run
    pool = multiprocessing.Pool(processes=n_processes)
    cluster_dicts = [
        pool.apply_async(populate_template_dict, args=(sequence_storage, f, sample_ids)) for f in fasta_files
    ]
    cluster_dicts = [result.get() for result in cluster_dicts]

    # Merge all of the dictionaries into one
    for d in cluster_dicts:
        for cluster, sequence in d.items():
            sequence_storage[cluster] += sequence
    return sequence_storage


def populate_template_dict(template_dict: dict, cluster_file: Path, sample_ids: list):
    with open(str(cluster_file), 'r') as f:
        cluster_dict = deepcopy(template_dict)
        seq_length = 0
        cluster_samples = []
        seq_lengths = []
        for line in f:
            if line.startswith(">"):
                # Add the sequence length to the ongoing tracking list
                if seq_length != 0:
                    seq_lengths.append(seq_length)
                seq_length = 0
                # Remove '>'
                header = line.replace(">", "").strip()
                # Add this header to cluster_samples to keep track of which samples we have sequence for
                cluster_samples.append(header)
            else:
                cluster_dict[header] += line.strip()
                seq_length += len(line)
        try:
            assert len(set(seq_lengths)) == 1
        except AssertionError:
            print(f"ERROR: Varying sequence lengths detected in {f.name}")
            print(seq_lengths)
        cluster_samples = set(cluster_samples)
        for s in sample_ids:
            if s not in cluster_samples:
                cluster_dict[s] += ('N' * seq_lengths[0])
        return cluster_dict
