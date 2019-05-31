from pathlib import Path

from dataclasses import dataclass


@dataclass
class SampleObject:
    """ Dataclass to store metadata for a sample """

    # Must be instantiated with these attributes
    sample_id: str
    fasta_path: Path

    # Updated later in the lifecycle
    prokka_object = None
    mmseqs_object = None
    cluster: Path = None
    sorted_faa: Path = None

    def __lt__(self, other):
        return self.sample_id < other.sample_id
