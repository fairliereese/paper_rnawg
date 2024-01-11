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
    combined_cds_path: str,
    source_gtf_path: str,
    cpat_cds_path: str,
    orfanage_cds_path: str,
    output_path: str,
) -> int:
    combined_cds = pd.read_csv(combined_cds_path, header=None, sep="\t")
    source_gtf = pd.read_csv(source_gtf_path, header=None, sep="\t")
    cpat_cds = pd.read_csv(cpat_cds_path, header=None, sep="\t")
    orfanage_cds = pd.read_csv(orfanage_cds_path, header=None, sep="\t")

    combined_cds["transcript_id"] = (
        combined_cds.iloc[:, GTF_GENE_INFO_IX]
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

    cpat_cds["transcript_id"] = (
        cpat_cds.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )

    orfanage_cds["transcript_id"] = (
        orfanage_cds.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    transcript_id_ix = 9
    mapping = source_gtf.iloc[:, [GTF_SOURCE_IX, transcript_id_ix]].drop_duplicates()
    mapping_cds = pd.concat([cpat_cds, orfanage_cds], axis=0).reset_index(drop=True)
    mapping_cds = (
        mapping_cds.loc[mapping_cds.iloc[:, GTF_TYPE_IX] == "CDS"]
        .iloc[:, [GTF_SOURCE_IX, transcript_id_ix]]
        .drop_duplicates()
    )
    combined_cds.loc[combined_cds.iloc[:, GTF_TYPE_IX] != "CDS", GTF_SOURCE_IX] = (
        combined_cds.loc[combined_cds.iloc[:, GTF_TYPE_IX] != "CDS"]
        .merge(mapping, on="transcript_id", how="left")["1_y"]
        .values
    )
    combined_cds.loc[combined_cds.iloc[:, GTF_TYPE_IX] == "CDS", GTF_SOURCE_IX] = (
        combined_cds.loc[combined_cds.iloc[:, GTF_TYPE_IX] == "CDS"]
        .merge(mapping_cds, on="transcript_id", how="left")["1_y"]
        .values
    )
    combined_cds = combined_cds.iloc[:, :-1]
    combined_cds.iloc[:, GTF_GENE_INFO_IX] = combined_cds.iloc[
        :, GTF_GENE_INFO_IX
    ].apply(lambda x: x.rsplit(" gene_name")[0])

    combined_cds.to_csv(
        f"{output_path}", sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE
    )
    return 0


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--combined_cds_path")
    parser.add_argument("--source_gtf_path")
    parser.add_argument("--cpat_cds_path")
    parser.add_argument("--orfanage_cds_path")
    parser.add_argument("--output_path")

    args = parser.parse_args()
    main(
        combined_cds_path=args.combined_cds_path,
        source_gtf_path=args.source_gtf_path,
        cpat_cds_path=args.cpat_cds_path,
        orfanage_cds_path=args.orfanage_cds_path,
        output_path=args.output_path,
    )
