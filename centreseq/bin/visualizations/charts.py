import numpy as np
import pandas as pd
import seaborn as sns

from random import sample
from pathlib import Path
from statistics import mean
from centreseq.bin.tree.summary import read_summary_report
from tqdm import tqdm


def generate_gene_count_dict(df: pd.DataFrame, repeats: int = 1):
    """
    Given n genomes
        1) How many clusters do they share (core)?
        2) How many total unique clusters do they have between them (pan)?

    - First, randomly pick n (start at 1) columns from df - remove rest
    - Next, answer questions 1) and 2)
    - Repeat with n+1 til we reach a maximum (this is equal to the total # of samples we have)

    Note that increasing the number of repeats increases runtime while tightening up the
    variance observed in the curves.

    """
    upper_range = df['n_members'].max()

    gene_count_dict = {x: {"core": [], "pan": []} for x in range(2, upper_range + 1)}
    for i in tqdm(range(2, upper_range + 1), desc="Calculating rarefaction"):
        # Repeat the process according to the number of repeats; end up with a mean value
        # TODO: maybe make box plots out of this data?
        for r in range(0, repeats):
            df_tmp = pick_random_columns(df, i).reset_index(drop=True)

            # Drop all completely empty rows then count total # rows remaining as the 'pan' genome
            df_pan = df_tmp.dropna(how='all')
            gene_count_dict[i]['pan'].append(df_pan.shape[0])

            # Drop any rows with any nan values present (i.e. not a core cluster among the subset)
            df_core = df_pan.dropna(how='any')
            gene_count_dict[i]['core'].append(df_core.shape[0])

        # Calculate mean value
        gene_count_dict[i]['core'] = mean(gene_count_dict[i]['core'])
        gene_count_dict[i]['pan'] = mean(gene_count_dict[i]['pan'])

    return gene_count_dict


def pick_random_columns(df: pd.DataFrame, n: int):
    """ Pick n random columns (excluding first 3 metadata columns), return filtered df """
    start = 4
    end = df.shape[1]
    rand_cols = sample(range(start, end), n)
    cols_to_select = rand_cols
    df_filtered = df.iloc[:, cols_to_select]
    return df_filtered


def generate_rarefaction_chart(summary_report: Path):
    outdir = summary_report.parent
    df_ = read_summary_report(summary_report)
    df_ = df_.replace(r'^\s*$', np.nan, regex=True)  # Replace empty strings with NaN values
    gene_count_dict_ = generate_gene_count_dict(df=df_, repeats=5)
    df_count = pd.DataFrame.from_dict(data=gene_count_dict_, orient='index')
    df_count.to_csv(outdir / 'rarefaction_curve.csv')
    sns.set(style="whitegrid", color_codes=True)
    ax = sns.scatterplot(data=df_count, linewidth=0, alpha=0.7)
    ax.set(xlabel="Number of genomes", ylabel="Number of genes")
    ax.figure.savefig(outdir / 'rarefaction_curve.png')
