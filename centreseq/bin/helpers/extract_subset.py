import sys
import time
from pathlib import Path

import pandas as pd

"""
Given an input text file of Sample IDs and a summary report, will return a filtered version of the summary report
of genes that belong exclusively in the input sample ID list and not the remainder
"""


def extract_subset(input_samples: Path, summary_report: Path, outpath: Path):
    if outpath is None:
        outpath = summary_report.parent / f'summary_report_subset_{time.strftime("%Y%m%d-%H%M%S")}.tsv'
    else:
        if outpath.suffix == "":
            outpath = outpath.with_suffix(".tsv")

    print(f"Started extract_subset.py with the following input values:")
    print(f"\tinput_samples\t: {input_samples}")
    print(f"\tsummary_report\t: {summary_report}")
    print(f"\toutpath\t\t\t: {outpath}")

    # Read in the user input
    df = pd.read_csv(summary_report, sep="\t")
    samples = read_input_samples(input_samples=input_samples)

    # Ensure that the user input valid samples
    validate_input_samples(samples=samples, summary_report_df=df)

    # Parse the DataFrame
    excluded_columns = ['cluster', 'cluster_representative', 'product', 'n_members']
    valid_columns = [x for x in list(df.keys()) if x not in excluded_columns]

    # Create a list of samples in the summary report not including the user's input samples
    other_samples = list(set(valid_columns) - set(samples))

    print(f"Filtering {summary_report.name}...")

    # Sort the dataframe so our desired rows are up top (this step isn't necessary, whatever)
    df = df.sort_values(by=other_samples, na_position='first')

    # Filter the dataframe, retaining only rows where all other_samples values are null
    df_filtered = df[(df[other_samples].isnull()).all(axis=1)]

    # Export df_filtered to csv
    df_filtered.to_csv(outpath, sep="\t", index=None)

    print(f"DONE! Output available at {outpath}")


def validate_input_samples(samples: list, summary_report_df: pd.DataFrame):
    excluded_columns = ['cluster', 'cluster_representative', 'product', 'n_members']
    valid_samples = [x for x in list(summary_report_df.keys()) if x not in excluded_columns]

    error_list = []
    for s in samples:
        if s not in valid_samples:
            error_list.append(s)
    if len(error_list) > 0:
        print(f"ERROR: The following input Sample IDs could not be found in the summary report:")
        print("\n".join(error_list))
        sys.exit()
    print("Successfully validated sample list")


def read_input_samples(input_samples: Path) -> list:
    samples = []
    with open(str(input_samples), 'r') as f:
        for line in f:
            samples.append(line.strip())
    return samples
