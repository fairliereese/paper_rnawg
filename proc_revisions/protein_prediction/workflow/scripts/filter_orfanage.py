#!/usr/bin/env python
# -*- coding: utf-8 -*-


from argparse import ArgumentParser

import numpy as np
import pandas as pd
from typeguard import typechecked

STOP_CODON_NT = 3

ZERO_ONE_INDEX_OFFSET = 1
GTF_TYPE_IX = 2
GTF_START_IX = 3
GTF_END_IX = 4
GTF_GENE_INFO_IX = 8


@typechecked
def main(
    input_file_path: str,
    output_path_to_be_predicted: str,
    output_path_filtered: str,
    minimum_orf_length: int,
) -> int:
    df = pd.read_csv(input_file_path, sep="\t", header=None)
    df["transcript_id"] = (
        df.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    transcript_id_index = len(df.columns) - 1
    cds = df.iloc[
        :, [GTF_TYPE_IX, GTF_START_IX, GTF_END_IX, transcript_id_index]
    ]
    cds.columns = ["type", "start", "stop", "transcript_id"]
    cds = cds.loc[cds["type"] == "CDS"]
    cds["width"] = cds["stop"] - cds["start"] + ZERO_ONE_INDEX_OFFSET

    cds["cumwidth"] = cds.groupby("transcript_id")["width"].transform(
        lambda x: pd.Series(x).cumsum()
    )
    cds_with_sufficient_length = set(
        np.unique(cds.loc[cds["cumwidth"] >= minimum_orf_length].transcript_id)
    )

    # Reload GTF to avoid cleanup.
    output_df = pd.read_csv(input_file_path, sep="\t", header=None)
    mask = np.array(
        [i in cds_with_sufficient_length for i in df.transcript_id.values]
    )
    output_df.loc[mask, :].to_csv(
        f"{output_path_filtered}",
        sep="\t",
        index=False,
        header=False,
        quotechar="'",
    )

    output_df = output_df.loc[np.logical_not(mask), :]
    output_df = output_df.loc[output_df.iloc[:, GTF_TYPE_IX] != "CDS"]
    output_df.to_csv(
        f"{output_path_to_be_predicted}",
        sep="\t",
        index=False,
        header=False,
        quotechar="'",
    )
    return 0


if __name__ == "__main__":

    parser = ArgumentParser()

    parser.add_argument("--orfanage_gtf_file_path")
    parser.add_argument("--output_path_to_be_predicted")
    parser.add_argument("--output_path_filtered")
    parser.add_argument("--minimum_orf_length")

    args = parser.parse_args()
    main(
        input_file_path=args.orfanage_gtf_file_path,
        output_path_to_be_predicted=args.output_path_to_be_predicted,
        output_path_filtered=args.output_path_filtered,
        minimum_orf_length=int(args.minimum_orf_length),
    )
