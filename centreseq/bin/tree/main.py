import logging
import multiprocessing
from pathlib import Path

import numpy as np
from tqdm import tqdm

from centreseq.bin.tree.cluster_data_structures import ClusterVariants, populate_cluster_object
from centreseq.bin.tree.concatenate_sequences import concatenate_sequence_directory
from centreseq.bin.tree.summary import generate_variant_summary_df, generate_variant_detail_df, write_tsv_from_df, \
    read_summary_report
from centreseq.bin.tree.wrappers import call_snp_sites, call_muscle

main_log = logging.getLogger('main_log')


def tree_pipeline(summary_report: Path, prokka_dir: Path, outdir: Path, n_cpu: int, percentile: float):
    """
    Main pipeline call to extract data, call SNPs
    :param summary_report: Report generated by core pipeline (e.g. summary_report.tsv)
    :param prokka_dir: Directory containing all Prokka results generated through core pipeline
    :param outdir: Output directory
    :param n_cpu: Number of CPUs to allocate to tasks that support multithreading
    :param percentile: Filter summary report to the top nth percentile
    """
    # Initial reading a filtering of the summary report
    df = read_summary_report(summary_report)
    original_length = len(df)

    # Filter out singleton values
    df = df[df['n_members'] > 1]

    # Retrieve member list
    member_list = [x for x in df.columns.values if x not in
                   ['cluster', 'cluster_representative', 'product', 'n_members']]

    # Filter to nth percentile of samples
    n_members_max = np.percentile(df['n_members'], percentile)
    main_log.debug(f"Filtering DataFrame to contain only rows where n_members >= {n_members_max}")
    df = df[df['n_members'] >= n_members_max]

    filtered_length = len(df)
    main_log.info(f"Extracting sequence information for core clusters described in {summary_report}")
    main_log.debug(f"Processing {filtered_length} clusters (filtered from {original_length} total clusters)")
    main_log.debug(f"{filtered_length}x{n_members_max} = {filtered_length * n_members_max} sequences will be extracted")

    # Extract sequence data for each cluster member of every cluster (this is slow)
    tqdm.pandas()
    cluster_objects = df.progress_apply(populate_cluster_object, axis=1, args=(prokka_dir,))

    # Write extracted sequences to a file for each cluster
    main_log.debug(f"Writing sequences to cluster files")
    aligned_loci_dir = outdir / 'aligned_loci'
    aligned_loci_dir.mkdir(exist_ok=True)
    for cluster_object in tqdm(cluster_objects):
        cluster_object.generate_cluster_fasta(outdir=aligned_loci_dir)

    # Align cluster multi-FASTA files
    main_log.debug(f"Aligning all cluster files with Muscle")
    p = multiprocessing.Pool(processes=n_cpu)
    cluster_fastas = [cluster_object.cluster_fasta for cluster_object in cluster_objects]
    for _ in tqdm(p.imap_unordered(func=call_muscle, iterable=cluster_fastas), total=len(cluster_fastas)):
        pass

    # Call variants with snp-sites, create a ClusterVariants object per variant
    """
    NOTE: With a large number of calls to snp-sites, there is a chance that the following error will occur:
        OSError [24]: Too many open files
    
    Not sure of a good permanent fix for this problem as I suspect that snp-sites is leaving some files open when 
    it shouldn't be.
    
    To circumvent this, I had to change a number of OS settings. Add the following lines to the following files:
        /etc/systemd/user.conf      DefaultLimitNOFILE=65535
        /etc/systemd/system.conf    DefaultLimitNOFILE=65535
        /etc/security/limits.conf   * soft nofile 65535
                                    * hard nofile 65535
                                    root soft nofile 65535
                                    root hard nofile 65535
    """
    main_log.debug(f"Calling variants for all aligned clusters")
    vcf_dir = outdir / 'variant_calls'
    vcf_dir.mkdir(exist_ok=True)
    variants_objects = []
    for cluster_object in tqdm(cluster_objects):
        vcf = call_snp_sites(cluster_object.cluster_fasta, vcf_dir)
        if vcf.exists():
            variants = ClusterVariants(parent_cluster=cluster_object, vcf_path=vcf)
            variants_objects.append(variants)
        else:
            variants = ClusterVariants(parent_cluster=cluster_object)
            variants_objects.append(variants)

    # Generate variant summary tables
    variant_summary_df = generate_variant_summary_df(cluster_variants_list=variants_objects)
    variant_detail_df = generate_variant_detail_df(cluster_variants_list=variants_objects)
    write_tsv_from_df(df=variant_summary_df, outpath=outdir / 'variants_summary.tsv')
    write_tsv_from_df(df=variant_detail_df, outpath=outdir / 'variants_detail.tsv')

    # Concatenate all core genes into a file to feed into RAXml
    main_log.debug(f"Concatenating core gene files")
    concatenated_sequence_dir = outdir / 'concatenated_sequences'
    concatenated_sequence_dir.mkdir(exist_ok=True)
    concatenate_sequence_directory(sample_ids=member_list, sequence_directory=aligned_loci_dir, n_processes=n_cpu,
                                   outdir=concatenated_sequence_dir)
