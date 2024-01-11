#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
from argparse import ArgumentParser

import pandas as pd
from typeguard import typechecked

GTF_SOURCE_IX = 1
GTF_TYPE_IX = 2
GTF_GENE_INFO_IX = 8


@typechecked
def main(
    cpat_gtf_path: str,
    source_gtf_path: str,
    output_path: str,
) -> int:
    cpat_gtf = pd.read_csv(cpat_gtf_path, header=None, sep="\t")
    source_gtf = pd.read_csv(source_gtf_path, header=None, sep="\t")

    cpat_gtf["transcript_id"] = (
        cpat_gtf.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )

    source_gtf["transcript_id"] = (
        source_gtf.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    transcript_id_ix = 9
    mapping = source_gtf.iloc[
        :, [GTF_SOURCE_IX, transcript_id_ix]
    ].drop_duplicates()
    cpat_gtf.loc[cpat_gtf.iloc[:, GTF_TYPE_IX] != "CDS", GTF_SOURCE_IX] = (
        cpat_gtf.loc[cpat_gtf.iloc[:, GTF_TYPE_IX] != "CDS"]
        .merge(mapping, on="transcript_id", how="left")["1_y"]
        .values
    )

    cpat_gtf.iloc[:, :-1].to_csv(
        f"{output_path}",
        sep="\t",
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE
    )
    return 0


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--cpat_gtf_path")
    parser.add_argument("--source_gtf_path")
    parser.add_argument("--output_path")

    args = parser.parse_args()
    main(
        cpat_gtf_path=args.cpat_gtf_path,
        source_gtf_path=args.source_gtf_path,
        output_path=args.output_path,
    )
