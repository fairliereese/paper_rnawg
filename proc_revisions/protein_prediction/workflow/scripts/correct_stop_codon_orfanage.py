#!/usr/bin/env python
# -*- coding: utf-8 -*-


from argparse import ArgumentParser

import numpy as np
import pandas as pd
from typeguard import typechecked

GTF_TYPE_IX = 2
GTF_START_IX = 3
GTF_END_IX = 4
GTF_STRAND_IX = 6
GTF_GENE_INFO_IX = 8


@typechecked
def main(
    orfanage_gtf_file_path: str,
    output_path: str,
) -> int:
    df = pd.read_csv(orfanage_gtf_file_path, sep="\t", header=None)
    df["transcript_id"] = (
        df.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )

    transcript_and_exon_mask = np.isin(
        df.index, df.loc[df.iloc[:, GTF_TYPE_IX] != "CDS"].index
    )
    cds = df.loc[df.iloc[:, GTF_TYPE_IX] == "CDS"].copy(deep=True)
    cds["width"] = cds.iloc[:, GTF_END_IX] - cds.iloc[:, GTF_START_IX] + 1
    sense_strand = (
        cds.loc[df.iloc[:, GTF_STRAND_IX] == "+"]
        .sort_values([0, GTF_END_IX], ascending=False)
        .copy(deep=True)
    )
    antisense_strand = (
        cds.loc[df.iloc[:, GTF_STRAND_IX] == "-"]
        .sort_values([0, GTF_START_IX], ascending=True)
        .copy(deep=True)
    )
    sense_strand["cumwidth"] = sense_strand.groupby("transcript_id")[
        "width"
    ].transform(lambda x: pd.Series(x).cumsum().shift(fill_value=0.0))

    antisense_strand["cumwidth"] = antisense_strand.groupby("transcript_id")[
        "width"
    ].transform(lambda x: pd.Series(x).cumsum().shift(fill_value=0.0))

    cds_sense_keep_mask = np.isin(
        df.index,
        sense_strand.loc[
            (sense_strand["cumwidth"] + sense_strand["width"]) > 3
        ].index,
    )
    cds_antisense_keep_mask = np.isin(
        df.index,
        antisense_strand.loc[
            (
                antisense_strand["cumwidth"].values
                + antisense_strand["width"].values
            )
            > 3
        ].index,
    )
    sense_indices = (
        sense_strand.loc[
            sense_strand["cumwidth"].values + sense_strand["width"].values > 3
        ]
        .groupby("transcript_id", as_index=False)[GTF_END_IX]
        .idxmax()
        .iloc[:, 1]
        .values
    )
    sense_offsets_index = (
        sense_strand.loc[
            (sense_strand["cumwidth"].values + sense_strand["width"].values)
            > 3
        ]
        .reset_index(drop=True)
        .groupby("transcript_id", as_index=False)[GTF_END_IX]
        .idxmax()
        .iloc[:, 1]
        .values
    )
    sense_offsets_value = (
        (
            sense_strand.loc[
                (
                    sense_strand["cumwidth"].values
                    + sense_strand["width"].values
                )
                > 3
            ]["cumwidth"]
            - 3
        ).values
    )[sense_offsets_index]

    antisense_indices = (
        antisense_strand.loc[
            (
                antisense_strand["cumwidth"].values
                + antisense_strand["width"].values
            )
            > 3
        ]
        .groupby("transcript_id", as_index=False)[GTF_START_IX]
        .idxmin()
        .iloc[:, 1]
        .values
    )
    antisense_offsets_index = (
        antisense_strand.loc[
            (
                antisense_strand["cumwidth"].values
                + antisense_strand["width"].values
            )
            > 3
        ]
        .reset_index(drop=True)
        .groupby("transcript_id", as_index=False)[GTF_START_IX]
        .idxmin()
        .iloc[:, 1]
        .values
    )
    antisense_offsets_value = (
        (
            antisense_strand.loc[
                (
                    antisense_strand["cumwidth"].values
                    + antisense_strand["width"].values
                )
                > 3
            ]["cumwidth"]
            - 3
        ).values
    )[antisense_offsets_index]

    df.iloc[sense_indices, GTF_END_IX] = (
        df.iloc[sense_indices, GTF_END_IX] + sense_offsets_value
    )
    df.iloc[antisense_indices, GTF_START_IX] = (
        df.iloc[antisense_indices, GTF_START_IX] - antisense_offsets_value
    )
    df = df.iloc[:, : GTF_GENE_INFO_IX + 1]
    df = df.loc[
        np.logical_or(
            np.logical_or(cds_sense_keep_mask, cds_antisense_keep_mask),
            transcript_and_exon_mask,
        )
    ]
    df.to_csv(
        f"{output_path}",
        sep="\t",
        index=False,
        header=False,
        quotechar="'",
    )
    return 0


if __name__ == "__main__":

    parser = ArgumentParser()

    parser.add_argument("--orfanage_gtf_file_path")
    parser.add_argument("--output_path")

    args = parser.parse_args()
    main(
        orfanage_gtf_file_path=args.orfanage_gtf_file_path,
        output_path=args.output_path,
    )
