from itertools import combinations
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def generate_pairwise_gene_match_report(summary_report: Path, outdir: Path) -> Path:
    """ Given a summary report, generates pairwise_gene_match_report.tsv and drops it into the outdir """
    df = generate_pairwise_gene_match_df(summary_report)
    report = export_pairwise_gene_match_report(df, outdir)
    return report


def export_pairwise_gene_match_report(df: pd.DataFrame, outdir: Path) -> Path:
    """ Exports a .tsv file of the dataframe created by generated_pairwise_gene_match_df() """
    outname = outdir / 'pairwise_gene_match_report.tsv'
    df.to_csv(outname, sep="\t", index=None)
    return outname


def generate_pairwise_gene_match_df(summary_report: Path) -> pd.DataFrame:
    """
    Given an existing summary report, compares every single potential pair of genomes and returns a DataFrame containing
    metrics on the number of shared genes between the pair.
    """

    # Improves reading performance with na_filter=False, na values in our case are just empty cells i.e. ""
    df = pd.read_csv(summary_report, sep="\t", na_filter=False)

    # Pull list of samples from the df
    keys_to_remove = ['cluster', 'cluster_representative', 'product', 'n_members']  # Column headers that we don't want
    samples = list(df.keys())
    samples = [sample for sample in samples if sample not in keys_to_remove]

    # Get all possible pairs of samples
    pairs = list(combinations(samples, 2))

    # Create pair metrics dictionary
    pair_dict = {}
    for pair in tqdm(pairs, desc="Pairwise comparisons"):
        # Subset the main df using only the columns for the current sample pair
        df_subset = df[[pair[0], pair[1]]].copy()

        # Create new "match" column where True is returned if both columns contain a gene, else False
        df_subset['match'] = (df_subset.iloc[:, 0] != '') & (df_subset.iloc[:, 1] != '')

        # Filter out rows where both samples are missing the gene
        df_subset = df_subset[df_subset.iloc[:, 0] != df_subset.iloc[:, 1]]

        # Get # of total genes between the two samples
        total_genes_among_pair = len(df_subset)

        # Extract value counts from the new 'match' column
        match_counting = df_subset['match'].value_counts().to_frame().transpose()

        # Get n_match
        if True not in match_counting:
            n_match = 0
        else:
            n_match = match_counting[True][0]

        # Get n_nomatch
        if False not in match_counting:
            n_nomatch = 0
        else:
            n_nomatch = match_counting[False][0]

        # Populate dictionary with the pair as the key
        pair_dict[pair] = {'n_match': n_match,
                           'n_nomatch': n_nomatch,
                           'match_percentage': (n_match / total_genes_among_pair) * 100}

    # Convert dict to pd.DataFrame with each row representing a pair
    match_df = pd.DataFrame.from_dict(pair_dict, orient='index')
    match_df = match_df.reset_index()
    match_df = match_df.rename(index=str, columns={"level_0": "sample_1",
                                                   "level_1": "sample_2"})
    return match_df
