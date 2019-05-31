import logging
from collections.__init__ import Counter
from pathlib import Path

import pandas as pd

from centreseq.bin.core.summary import cluster_tsv_to_df, get_cluster_member_dict

log = logging.getLogger('main_log')


def filter_core_cluster_tsv(cluster_tsv: Path, outdir: Path):
    """
    Filters/repairs the core cluster .tsv output file to ensure there are no cases where a single sample is represented
    in a cluster more than once.

    :param cluster_tsv: path to cluster .tsv file
    :param outdir: path desired output directory
    :return: path to filtered/repaired cluster .tsv
    """

    # Read TSV
    df = cluster_tsv_to_df(cluster_tsv=cluster_tsv)

    # TODO: Factor some of the logic below into clean functions
    # Get the cluster member dict
    cluster_member_dict = get_cluster_member_dict(df)

    # Determine which clusters contain duplicate samples, create a dictionary to store them
    to_repair_dict = {}
    for cluster_rep, members in cluster_member_dict.items():
        sample_id_list = []
        for m in members:
            sample_id = m.rsplit("_", 1)[0]
            sample_id_list.append(sample_id)
        counter_dict = Counter(sample_id_list)
        problematic_sample_id_dict = {x: counter_dict[x] for x in counter_dict if counter_dict[x] > 1}
        if len(problematic_sample_id_dict) >= 1:
            log.debug(f"Warning: Detected problematic samples: {problematic_sample_id_dict}")
            to_repair_dict[cluster_rep] = problematic_sample_id_dict
        else:
            continue

    # Write a file containing all clusters/members flagged as problematic
    debug_dict = {}
    for cluster_rep, members in cluster_member_dict.items():
        sample_id_list = []
        for m in members:
            sample_id = m.rsplit("_", 1)[0]
            sample_id_list.append(sample_id)
        counter_dict = Counter(sample_id_list)
        problematic_sample_id_dict = {x: counter_dict[x] for x in counter_dict if counter_dict[x] > 1}
        if len(problematic_sample_id_dict) >= 1:
            debug_list = []
            for s in list(problematic_sample_id_dict.keys()):
                for m in members:
                    if s in m:
                        debug_list.append(m)
            debug_dict[cluster_rep] = debug_list

    debug_file = Path(outdir / 'master_genome_DB.cluster.debug.tsv')
    with open(str(debug_file), 'w') as f:
        for cluster_rep, members in debug_dict.items():
            for m in members:
                f.write(f"{cluster_rep}\t{m}\n")

    # Repair summary tsv by removing rows til only 1 gene per sample is represented within a cluster
    # Determine index of rows to filter out
    row_indexes_to_remove = []
    for i, row in enumerate(df.itertuples()):
        cluster_rep = row.cluster_rep
        cluster_member = row.cluster_member
        if cluster_rep in to_repair_dict.keys():
            # Member sample_IDs that need to be reduced
            members_dict = to_repair_dict[cluster_rep]
            members_list = list(members_dict.keys())

            # Detect if this is a row we want to remove
            for m in members_list:
                if m in cluster_member:
                    # Check counter dict to see if there is > 1 member for cluster
                    if to_repair_dict[cluster_rep][m] > 1:
                        log.debug(f"Repairing {cluster_rep}")
                        row_indexes_to_remove.append(i)
                        # Decrement the counter dict
                        to_repair_dict[cluster_rep][m] -= 1
        else:
            continue

    # Create filtered version of df with duplicate rows removed
    filtered_df = df.drop(df.index[row_indexes_to_remove])

    # Export to TSV file
    outpath = outdir / cluster_tsv.with_suffix(".filtered.tsv").name
    filtered_df.to_csv(outpath, sep="\t", header=None, index=None)
    get_difference_between_cluster_tsvs(cluster_tsv=cluster_tsv, filtered_cluster_tsv=outpath)
    return outpath


def get_difference_between_cluster_tsvs(cluster_tsv: Path, filtered_cluster_tsv: Path) -> Path:
    """ Creates a dataframe that only contains rows that differ between the two cluster dataframes """
    df1 = cluster_tsv_to_df(cluster_tsv)
    df2 = cluster_tsv_to_df(filtered_cluster_tsv)
    df = pd.concat([df1, df2]).drop_duplicates(keep=False)
    outpath = cluster_tsv.parent / cluster_tsv.with_suffix(".removed_rows.tsv").name
    df.to_csv(outpath, sep="\t", header=None, index=None)
    return outpath
