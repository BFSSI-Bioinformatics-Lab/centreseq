from pathlib import Path

import pandas as pd


def write_tsv_from_df(df: pd.DataFrame, outpath: Path) -> Path:
    """ Wrapper for pd.DataFrame.to_csv with preferred defaults """
    df.to_csv(outpath, sep="\t", index=None)
    return outpath


def read_summary_report(summary_report: Path) -> pd.DataFrame:
    """ Wrapper for pd.DataFrame.read_csv with preferred defaults """
    return pd.read_csv(summary_report, sep="\t", na_filter=False)
