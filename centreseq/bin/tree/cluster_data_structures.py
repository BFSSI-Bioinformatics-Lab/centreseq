import logging
from pathlib import Path

import pandas as pd
from cyvcf2.cyvcf2 import VCF
from dataclasses import dataclass

main_log = logging.getLogger('main_log')


@dataclass
class Cluster:
    """ Dataclass to store data for a core genome cluster """
    cluster_id: str
    cluster_representative: str
    product: str
    n_members: int
    member_list: list
    prokka_dir: Path
    vcf_path: Path = None
    cluster_fasta: Path = None

    def __post_init__(self):
        # Produce ClusterMember dataclass objects
        self.member_objects = self.generate_member_objects()

        # Retrieve sequence length by iterating through each member object and finding an entry that is not None
        self.cluster_sequence_length = self.get_cluster_sequence_length()

        # Iterate through member objects and replace None sequences with 'N' nucleotide padding
        self.pad_empty_member_sequences()

    def pad_empty_member_sequences(self):
        for o in self.member_objects:
            if o.sequence_length is None:
                o.populate_empty_sequence(sequence_length=self.cluster_sequence_length)

    def get_cluster_sequence_length(self):
        for o in self.member_objects:
            if o.sequence_length is not None:
                return o.sequence_length

    def generate_member_objects(self):
        """ Creates a ClusterMember instance for each entry in self.member_list """
        member_objects = []
        for member in self.member_list:
            # With a properly annotated Prokka sequence, member_id and locus_id will be extracted correctly
            member_id = member.rsplit("_", 1)[0]
            locus_id = member.rsplit("_", 1)[1]
            # Check if locus_id has been annotated with 'EMPTY-MEMBER' which indicates member has no match
            if locus_id == "EMPTY-MEMBER":
                member_object = ClusterMember(parent_cluster=self,
                                              locus_id=locus_id,
                                              member_id=member_id)
            else:
                member_object = ClusterMember(parent_cluster=self,
                                              prokka_faa=self.prokka_dir / member_id / f"{member_id}.faa",
                                              prokka_ffn=self.prokka_dir / member_id / f"{member_id}.ffn",
                                              locus_id=locus_id,
                                              member_id=member_id)
            member_objects.append(member_object)
        return member_objects

    def generate_cluster_fasta(self, outdir: Path):
        """ Writes a multi-FASTA file containing sequences for each cluster member """
        cluster_fasta = outdir / f"{self.cluster_id}.ffn"
        if cluster_fasta.exists():
            cluster_fasta.unlink()
        with open(str(cluster_fasta), "a+") as f:
            [f.write(f">{member.member_id}\n{member.ffn_sequence}\n") for member in self.member_objects]
        self.cluster_fasta = cluster_fasta


@dataclass
class ClusterMember:
    """ Dataclass to store metadata on each member of a Cluster, including the full nucleotide sequence """
    parent_cluster: Cluster
    member_id: str

    locus_id: str = None
    prokka_faa: Path = None  # Not sure if we need to do this or not, will ignore for now
    faa_sequence: str = None
    prokka_ffn: Path = None
    ffn_sequence: str = None
    sequence_length: int = None

    def __post_init__(self):
        # Retrieve nucleotide sequence from Prokka ffn file for input contig
        if self.prokka_ffn is not None:
            self.ffn_sequence = extract_sequence(self.prokka_ffn, target_contig=f"{self.member_id}_{self.locus_id}")
            # This value should be the same for all members of the Cluster
            self.sequence_length = len(self.ffn_sequence)

    def populate_empty_sequence(self, sequence_length):
        """ Method to populate the ffn_sequence with'N' nucleotides in lieu of actual sequence """
        self.ffn_sequence = "N" * sequence_length
        self.sequence_length = sequence_length


@dataclass
class ClusterVariants:
    """
    Dataclass to store variant output from a parsed VCF of a Cluster
    Variants are generated with snp-sites and parsed with cyvcf2

    Note on self.variant_list:
        Contains {Variant} objects - probably easiest to extract a per-sample genotype via the .genotypes attribute.
        Use the .ALT and .REF to determine what the index from .genotypes indicates
    """
    parent_cluster: Cluster
    vcf_path: Path = None

    def __post_init__(self):
        # Parse the VCF file and extract relevant information to store in self.variant_sample_ids and self.variant_list
        self.variant_list = []
        # Only try to parse VCF if the file was generated by snp-sites
        if self.vcf_path is not None:
            vcf = VCF(str(self.vcf_path))
            self.variant_sample_ids = vcf.samples
            for v in vcf:
                self.variant_list.append(v)
        else:
            self.variant_sample_ids = []

    @property
    def num_variants(self):
        return len(self.variant_list)

    @property
    def num_samples(self):
        return len(self.variant_sample_ids)

    def __len__(self):
        return len(self.variant_list)


def populate_cluster_object(row: pd.Series, prokka_dir: Path):
    """
    Creates a Cluster object from a row extracted from the summary_report.tsv file produced by the core pipeline
    :param row:
    :param prokka_dir: Path to the
    :return:
    """
    member_cols = [x for x in row.keys() if x not in ['cluster', 'cluster_representative', 'product', 'n_members']]

    """
    Iterate over columns containing member data and extract the member. 
    If the value for the column is missing (i.e. sample is not represented in the cluster), add the name of sample 
    with no loci ID
    """
    member_id_list = []
    for member_id in member_cols:
        if row[member_id] != "":
            member_id_list.append(row[member_id])
        else:
            # Add the "EMPTY-MEMBER" locus ID to detect later
            member_id_list.append(member_id + "_EMPTY-MEMBER")

    cluster_object = Cluster(cluster_id=row['cluster'], cluster_representative=row['cluster_representative'],
                             product=row['product'], n_members=row['n_members'], member_list=member_id_list,
                             prokka_dir=prokka_dir)
    return cluster_object


def extract_sequence(fasta: Path, target_contig: str) -> str:
    """
    Given an input FASTA file and a target contig, will extract the nucleotide/amino acid sequence
    and return as a string. Will break out after capturing 1 matching contig.
    """
    sequence = ""
    write_flag = False
    write_counter = 0
    with open(str(fasta), 'r') as infile:
        for line in infile:
            if write_counter > 1:
                break
            if line[0] == ">":  # This is marginally faster than .startswith()
                if target_contig in line:
                    write_counter += 1
                    write_flag = True
                else:
                    write_flag = False
            elif write_flag:
                sequence += line
            else:
                continue
    return sequence.replace("\n", "")
