from collections.__init__ import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from centreseq.bin.core.accessories import get_fasta_headers


def cluster_tsv_to_df(cluster_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(cluster_tsv, sep="\t", header=None, names=("cluster_rep", "cluster_member"), low_memory=False)
    return df


def get_member_counts(df: pd.DataFrame) -> dict:
    # Get a complete list of all cluster_rep values including duplicates, do total counts
    raw_cluster_rep_list = list(df['cluster_rep'])
    member_counts = Counter(raw_cluster_rep_list)
    return member_counts


def get_unique_cluster_rep_list(df: pd.DataFrame) -> list:
    # Get a unique list of cluster_reps
    unique_cluster_list = df['cluster_rep'].unique()
    return unique_cluster_list


def init_summary_df(unique_cluster_list) -> pd.DataFrame:
    # Create a cluster ID for each cluster_rep
    cluster_id_list = [f"CLUSTER{x:06}" for x in range(len(unique_cluster_list))]

    # Assign IDs to cluster_reps
    cluster_dict = dict(zip(cluster_id_list, unique_cluster_list))

    # Create new df
    df = pd.DataFrame.from_dict(cluster_dict, orient="index").reset_index()
    df.columns = ['cluster', 'cluster_representative']

    # Add empty product column - will be populated later
    df['product'] = ""
    return df


def retrieve_member_count(value, member_count_dict) -> dict:
    return member_count_dict[value]


def get_cluster_member_dict(initial_df: pd.DataFrame):
    """ Keys are cluster_rep values and values are a list of members for the cluster """
    data = initial_df.values
    cluster_member_dict = {d[0]: list() for d in data}
    for d in data:
        cluster_member_dict[d[0]].append(d[1])
    return cluster_member_dict


def get_sample_name(val) -> str:
    """ Returns the sample name from a prokka annotation e.g. BMH-2018-000398_00836 will return BMH-2018-000398 """
    return val.rsplit("_", 1)[0]


def get_samples_present(unique_cluster_list) -> set:
    samples_present = set([get_sample_name(val) for val in unique_cluster_list])
    return samples_present


def clust_rep_populate_with_apply(row, cluster_member_dict: dict):
    clust_rep = row['cluster_representative']
    clust_members = cluster_member_dict[clust_rep]
    for clust_member in clust_members:
        sample_name = get_sample_name(clust_member)
        row[sample_name] = clust_member
    return row


def add_members_to_df(initial_df, summary_df, sample_list) -> pd.DataFrame:
    cluster_member_dict = get_cluster_member_dict(initial_df=initial_df)

    # Prepopulate df with empty columns for each sample
    for sample in sample_list:
        summary_df[sample] = None

    # Populate df with cluster members
    summary_df = summary_df.apply(func=clust_rep_populate_with_apply, cluster_member_dict=cluster_member_dict, axis=1)
    return summary_df


def get_rep_seq_product_dict(fasta_header_list: list) -> dict:
    rep_seq_product_dict = {}
    for header in fasta_header_list:
        representative_seq = header.split(" ", 1)[0].replace(">", "")
        product = header.split(" ", 1)[1]
        rep_seq_product_dict[representative_seq] = product
    return rep_seq_product_dict


def add_product_column(summary_df: pd.DataFrame, core_genome: Path) -> pd.DataFrame:
    fasta_header_list = get_fasta_headers(core_genome)
    rep_seq_product_dict = get_rep_seq_product_dict(fasta_header_list)
    for index, row in summary_df.iterrows():
        clust_rep = row['cluster_representative']
        summary_df.at[index, 'product'] = rep_seq_product_dict[clust_rep]
    return summary_df


def export_to_tsv(outdir: Path, df: pd.DataFrame, filename: str) -> Path:
    outname = outdir / filename
    df.to_csv(outname, sep="\t", header=True, index=False)
    return outname


def generate_roary_gene_count_dict(summary_report_tsv: Path, n_samples: int) -> dict:
    """
    Emulates Roary default summary_statistics.txt
    """

    # [min, max]
    threshold_dict = {
        'core_genes': [int(0.99 * n_samples), n_samples],
        'soft_core_genes': [int(0.95 * n_samples), int(0.99 * n_samples)],
        'shell_genes': [int(0.15 * n_samples), int(0.95 * n_samples)],
        'cloud_genes': [0, int(0.15 * n_samples)],
        'total_genes': [0, n_samples]
    }

    df = pd.read_csv(summary_report_tsv, sep='\t', low_memory=False)
    core_genes = df[(df.n_members >= threshold_dict['core_genes'][0]) &
                    (df.n_members <= threshold_dict['core_genes'][1])]
    n_core_genes = len(core_genes)

    soft_core_genes = df[(df.n_members >= threshold_dict['soft_core_genes'][0]) &
                         (df.n_members < threshold_dict['soft_core_genes'][1])]
    n_soft_core_genes = len(soft_core_genes)

    shell_genes = df[(df.n_members >= threshold_dict['shell_genes'][0]) &
                     (df.n_members < threshold_dict['shell_genes'][1])]
    n_shell_genes = len(shell_genes)

    cloud_genes = df[(df.n_members >= threshold_dict['cloud_genes'][0]) &
                     (df.n_members < threshold_dict['cloud_genes'][1])]
    n_cloud_genes = len(cloud_genes)

    n_total_genes = len(
        df[(df.n_members >= threshold_dict['total_genes'][0]) &
           (df.n_members <= threshold_dict['total_genes'][1])])

    roary_gene_df_dict = {
        'core_genes': core_genes,
        'soft_core_genes': soft_core_genes,
        'shell_genes': shell_genes,
        'cloud_genes': cloud_genes,
    }

    roary_gene_count_dict = {
        'core_genes': n_core_genes,
        'soft_core_genes': n_soft_core_genes,
        'shell_genes': n_shell_genes,
        'cloud_genes': n_cloud_genes,
        'total_genes': n_total_genes
    }
    return roary_gene_count_dict


def generate_core_gene_count_dict(summary_report_tsv: Path, n_samples: int) -> dict:
    threshold_dict = {
        'hard_core': 1.0 * n_samples,
        'medium_core': int(0.90 * n_samples),
        'soft_core': int(0.50 * n_samples),
    }

    df = pd.read_csv(summary_report_tsv, sep='\t', low_memory=False)

    # Chunk df into hard, medium, soft
    n_hard_core_genes = len(df[df.n_members >= threshold_dict['hard_core']])
    n_medium_core_genes = len(df[df.n_members >= threshold_dict['medium_core']])
    n_soft_core_genes = len(df[df.n_members >= threshold_dict['soft_core']])

    core_gene_count_dict = {
        'n_hard_core_genes': n_hard_core_genes,
        'n_medium_core_genes': n_medium_core_genes,
        'n_soft_core_genes': n_soft_core_genes
    }

    return core_gene_count_dict


def generate_roary_gene_count_report(roary_gene_count_dict: dict, outdir: Path, coverage_length: float,
                                     min_seq_id: float) -> Path:
    """
    Equivalent of summary_statistics.txt from a Roary run
    """
    roary_gene_count_report = outdir / 'roary_gene_count_report.txt'
    if roary_gene_count_report.exists():
        roary_gene_count_report.unlink()

    with open(str(roary_gene_count_report), 'w') as f:
        f.write(f"Core genes\t(99% <= strains <= 100%)\t{roary_gene_count_dict['core_genes']}\n")
        f.write(f"Soft core genes\t(95% <= strains < 99%)\t{roary_gene_count_dict['soft_core_genes']}\n")
        f.write(f"Shell genes\t(15% <= strains < 95%)\t{roary_gene_count_dict['shell_genes']}\n")
        f.write(f"Cloud genes\t(0% <= strains < 15%)\t{roary_gene_count_dict['cloud_genes']}\n")
        f.write(f"Total genes\t(0% <= strains <= 100%)\t{roary_gene_count_dict['total_genes']}")
        f.write(
            f"\n\nResults were generated with:\n\tmin_seq_id\t= {min_seq_id}\n\tcoverage_length\t= {coverage_length}")

    return roary_gene_count_report


def generate_core_gene_count_report(core_gene_count_dict: dict, outdir: Path, coverage_length: float,
                                    min_seq_id: float) -> Path:
    core_gene_count_report = outdir / 'core_gene_count_report.txt'
    if core_gene_count_report.exists():
        core_gene_count_report.unlink()

    with open(str(core_gene_count_report), 'w') as f:
        f.write(f"# genes shared among 100% of samples:\t{core_gene_count_dict['n_hard_core_genes']}\n")
        f.write(f"# genes shared in >=90% of samples:\t{core_gene_count_dict['n_medium_core_genes']}\n")
        f.write(f"# genes shared in >=50% of samples:\t{core_gene_count_dict['n_soft_core_genes']}")
        f.write(
            f"\n\nResults were generated with:\n\tmin_seq_id\t= {min_seq_id}\n\tcoverage_length\t= {coverage_length}")

    return core_gene_count_report


def remove_singletons_summary_report(df: pd.DataFrame) -> pd.DataFrame:
    """
    Removes singleton rows from summary report
    """
    df = df.drop(df[df['n_members'] <= 1].index)
    return df


def generate_summary_report_dataframe(cluster_tsv: Path, core_genome: Path):
    df_ = cluster_tsv_to_df(cluster_tsv)
    member_count_dict = get_member_counts(df_)
    unique_cluster_list = get_unique_cluster_rep_list(df_)
    df = init_summary_df(unique_cluster_list)
    df = add_product_column(summary_df=df, core_genome=core_genome)
    df['n_members'] = np.vectorize(retrieve_member_count)(df['cluster_representative'], member_count_dict)
    samples_list = get_samples_present(unique_cluster_list)
    df = add_members_to_df(initial_df=df_, summary_df=df, sample_list=samples_list)
    df = df.sort_values(by=['n_members'], ascending=False)
    return df


def generate_summary_report(cluster_tsv: Path, core_genome: Path, outdir: Path) -> tuple:
    df = generate_summary_report_dataframe(cluster_tsv=cluster_tsv, core_genome=core_genome)
    filtered_df = remove_singletons_summary_report(df=df)
    report_tsv_path = export_to_tsv(outdir=outdir, df=df, filename="summary_report.tsv")
    report_tsv_path_filtered = export_to_tsv(outdir, df=filtered_df, filename="summary_report_singletons_removed.tsv")
    return report_tsv_path, report_tsv_path_filtered
