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
GTF_STRAND_IX = 6
GTF_GENE_INFO_IX = 8


@typechecked
def main(
    transcript_exons_path: str,
    cds_only_path: str,
    output_prefix: str,
) -> int:
    transcript = pd.read_csv(
        f"{transcript_exons_path}",
        sep="\t",
        header=None,
    )
    transcript = transcript.loc[transcript.iloc[:, GTF_TYPE_IX] != "transcript"]
    cds = pd.read_csv(
        f"{cds_only_path}",
        sep="\t",
        header=None,
    )
    cds = cds.loc[cds.iloc[:, GTF_TYPE_IX] != "transcript"]
    transcript["transcript_id"] = (
        transcript.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
    )

    cds["transcript_id"] = (
        cds.iloc[:, GTF_GENE_INFO_IX].str.rsplit(";").str[0].str.rsplit('"').str[1]
    )
    transcript["width"] = (
        transcript.iloc[:, GTF_END_IX]
        - transcript.iloc[:, GTF_START_IX]
        + ZERO_ONE_INDEX_OFFSET
    )
    cds["width"] = (
        cds.iloc[:, GTF_END_IX] - cds.iloc[:, GTF_START_IX] + ZERO_ONE_INDEX_OFFSET
    )
    transcript["cumwidth"] = transcript.groupby("transcript_id")["width"].transform(
        lambda x: pd.Series(x).cumsum().shift(fill_value=0.0) + 1
    )
    cds["cumwidth"] = cds.groupby("transcript_id")["width"].transform(
        lambda x: pd.Series(x).cumsum() + 2
    )

    cds_merge = (
        cds.groupby("transcript_id")
        .head(1)
        .iloc[:, :-1]
        .sort_values("transcript_id")
        .reset_index(drop=True)
    )
    cds_merge["cumwidth"] = (
        cds.sort_values("cumwidth", ascending=False)
        .groupby("transcript_id")
        .head(1)
        .sort_values("transcript_id")["cumwidth"]
        .reset_index(drop=True)
    )

    cds_merge = cds_merge.iloc[:, [GTF_START_IX, GTF_STRAND_IX, 9, 11]]
    cds_merge.columns = ["cds_start", "strand", "transcript_id", "cumwidth"]
    transcript.transcript_id = transcript.transcript_id.astype(str)
    cds_merge.transcript_id = cds_merge.transcript_id.astype(str)

    combined = pd.merge(
        transcript,
        cds_merge,
        left_on="transcript_id",
        right_on="transcript_id",
    )
    combined["cds_start_transcript_space"] = (
        combined.iloc[:, GTF_START_IX] - combined["cds_start"]
    )
    range_finder = combined.loc[combined["cds_start_transcript_space"] <= 0]

    annotation = pd.DataFrame(
        {
            "transcript_id": range_finder.sort_values(
                ["transcript_id", "cds_start_transcript_space"],
                ascending=False,
            )
            .groupby("transcript_id")
            .head(1)["transcript_id"]
            .values,
            "orf_start": range_finder.sort_values(
                ["transcript_id", "cds_start_transcript_space"],
                ascending=False,
            )
            .groupby("transcript_id")
            .head(1)["cumwidth_x"]
            .values
            - range_finder.sort_values(
                ["transcript_id", "cds_start_transcript_space"],
                ascending=False,
            )
            .groupby("transcript_id")
            .head(1)["cds_start_transcript_space"]
            .values,
            "orf_end": range_finder.sort_values(
                ["transcript_id", "cds_start_transcript_space"],
                ascending=False,
            )
            .groupby("transcript_id")
            .head(1)["cumwidth_x"]
            .values
            - range_finder.sort_values(
                ["transcript_id", "cds_start_transcript_space"],
                ascending=False,
            )
            .groupby("transcript_id")
            .head(1)["cds_start_transcript_space"]
            .values
            + range_finder.sort_values(
                ["transcript_id", "cds_start_transcript_space"],
                ascending=False,
            )
            .groupby("transcript_id")
            .head(1)["cumwidth_y"]
            .values,
            "strand": range_finder.sort_values(
                ["transcript_id", "cds_start_transcript_space"],
                ascending=False,
            )
            .groupby("transcript_id")
            .head(1)["strand"]
            .values,
        }
    )

    annotation["orf_frame"] = "/"
    annotation["orf_len"] = (
        annotation["orf_end"] - annotation["orf_start"] + ZERO_ONE_INDEX_OFFSET
    )

    transcript = pd.read_csv(
        f"{output_prefix}.transcript_exons_only.gtf", sep="\t", header=None
    )
    transcript = transcript.loc[transcript.iloc[:, GTF_TYPE_IX] != "transcript"]
    transcript["transcript_id"] = (
        transcript.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
    )
    set_mask = set(
        np.intersect1d(transcript["transcript_id"], annotation["transcript_id"])
    )
    transcript = transcript.loc[
        np.array([i in set_mask for i in transcript.transcript_id]), :
    ]
    transcript["width"] = (
        transcript.iloc[:, GTF_END_IX]
        - transcript.iloc[:, GTF_START_IX]
        + ZERO_ONE_INDEX_OFFSET
    )
    transcript["cumwidth"] = transcript.groupby("transcript_id")["width"].transform(
        lambda x: pd.Series(x).cumsum()
    )

    annotation["len"] = (
        transcript.sort_values(["transcript_id", "cumwidth"], ascending=False)
        .groupby("transcript_id")
        .head(1)["cumwidth"]
        .values
    )

    annotation_positive_strand = annotation.loc[annotation["strand"] == "+"].copy(
        deep=True
    )
    annotation_negative_strand = annotation.loc[annotation["strand"] == "-"].copy(
        deep=True
    )
    annotation_negative_strand["orf_end"] = (
        annotation_negative_strand["len"]
        - annotation_negative_strand["orf_start"]
        + STOP_CODON_NT
        + ZERO_ONE_INDEX_OFFSET
    )
    annotation_negative_strand["orf_start"] = (
        annotation_negative_strand["orf_end"]
        - annotation_negative_strand["orf_len"]
        + ZERO_ONE_INDEX_OFFSET
    )

    annotation = pd.concat(
        [annotation_positive_strand, annotation_negative_strand],
        axis=0,
        ignore_index=True,
    )

    annotation[
        [
            "transcript_id",
            "len",
            "orf_frame",
            "orf_start",
            "orf_end",
            "orf_len",
        ]
    ].to_csv(f"{output_prefix}_best_orf.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--transcript_exons_path")
    parser.add_argument("--cds_only_path")
    parser.add_argument("--output_prefix")

    args = parser.parse_args()
    main(
        transcript_exons_path=args.transcript_exons_path,
        cds_only_path=args.cds_only_path,
        output_prefix=args.output_prefix,
    )
