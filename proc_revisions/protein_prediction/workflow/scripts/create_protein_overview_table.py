#!/usr/bin/env python
# -*- coding: utf-8 -*-


from argparse import ArgumentParser

import numpy as np
import pandas as pd
from Bio import SeqIO
from typeguard import typechecked

CODON_LENGTH = 3
STOP_CODON_NT = 3

GTF_SOURCE_IX = 1
GTF_TYPE_IX = 2
GTF_GENE_INFO_IX = 8

GTF_START_IX = 3
GTF_END_IX = 4
GTF_STRAND_IX = 6


@typechecked
def main(
    best_orf_path: str,
    sqanti_protein_path: str,
    orf_completeness_path: str,
    output_name: str,
    gtf_original_path: str,
    gtf_predicted_path: str,
    protein_fasta_path: str,
    blastp_path: str,
) -> int:
    orfs = pd.read_csv(
        f"{best_orf_path}",
        sep="\t",
    ).sort_values("transcript_id")
    orf_completeness = pd.read_csv(
        f"{orf_completeness_path}",
        sep="\t",
    ).sort_values("transcript_id")

    sqanti_protein = pd.read_csv(
        f"{sqanti_protein_path}",
        sep="\t",
    ).sort_values("pb")

    common_transcript_ids = sqanti_protein.pb.values

    orfs = (
        orfs.loc[np.isin(orfs.transcript_id, common_transcript_ids)]
        .copy(deep=True)
        .sort_values("transcript_id")
    )
    orf_completeness = (
        orf_completeness.loc[
            np.isin(orf_completeness.transcript_id, common_transcript_ids)
        ]
        .copy(deep=True)
        .sort_values("transcript_id")
    )
    sqanti_protein = sqanti_protein.copy(deep=True).sort_values("pb")
    gtf = pd.read_csv(f"{gtf_original_path}", sep="\t", header=None)
    gtf = gtf.loc[gtf.iloc[:, GTF_TYPE_IX] == "transcript"]
    gtf["transcript_id"] = (
        gtf.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[10]
        .str.rsplit('"')
        .str[1]
        .values
    )
    gtf = gtf.loc[np.isin(gtf.transcript_id.values, common_transcript_ids)]
    gtf["gene_id"] = (
        gtf.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    gtf["gene_name"] = (
        gtf.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[12]
        .str.rsplit('"')
        .str[1]
        .values
    )
    gtf = gtf.sort_values("transcript_id")
    gtf_predicted = pd.read_csv(f"{gtf_predicted_path}", sep="\t", header=None)
    cds_source = gtf_predicted.copy(deep=True)
    cds_source = cds_source.loc[(cds_source.iloc[:, GTF_TYPE_IX] == "CDS")]
    cds_source["transcript_id"] = (
        cds_source.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    cds_source = cds_source.iloc[:, [GTF_SOURCE_IX, -1]].drop_duplicates()
    cds_source.columns = ["source", "transcript_id"]
    cds_source = cds_source.sort_values("transcript_id")
    cds_positive_start = gtf_predicted.copy(deep=True)
    cds_positive_start["transcript_id"] = (
        cds_positive_start.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    cds_positive_start = cds_positive_start.loc[
        (cds_positive_start.iloc[:, GTF_TYPE_IX] == "CDS")
        & (cds_positive_start.iloc[:, GTF_STRAND_IX] == "+")
    ]
    cds_positive_start = (
        cds_positive_start.groupby("transcript_id", group_keys=True)
        .min()
        .iloc[:, GTF_START_IX]
        .reset_index()
        .sort_values("transcript_id")
    )
    cds_positive_start.columns = [cds_positive_start.columns[0], "coord"]

    cds_positive_end = gtf_predicted.copy(deep=True)
    cds_positive_end["transcript_id"] = (
        cds_positive_end.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    cds_positive_end = cds_positive_end.loc[
        (cds_positive_end.iloc[:, GTF_TYPE_IX] == "CDS")
        & (cds_positive_end.iloc[:, GTF_STRAND_IX] == "+")
    ]
    cds_positive_end = (
        cds_positive_end.groupby("transcript_id", group_keys=True)
        .max()
        .iloc[:, GTF_END_IX]
        .reset_index()
        .sort_values("transcript_id")
    )
    cds_positive_end.columns = [cds_positive_start.columns[0], "coord"]

    cds_negative_start = gtf_predicted.copy(deep=True)
    cds_negative_start["transcript_id"] = (
        cds_negative_start.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    cds_negative_start = cds_negative_start.loc[
        (cds_negative_start.iloc[:, GTF_TYPE_IX] == "CDS")
        & (cds_negative_start.iloc[:, GTF_STRAND_IX] == "-")
    ]
    cds_negative_start = (
        cds_negative_start.groupby("transcript_id", group_keys=True)
        .min()
        .iloc[:, GTF_START_IX]
        .reset_index()
        .sort_values("transcript_id")
    )
    cds_negative_start.columns = [cds_positive_start.columns[0], "coord"]

    cds_negative_end = gtf_predicted.copy(deep=True)

    cds_negative_end["transcript_id"] = (
        cds_negative_end.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    cds_negative_end = cds_negative_end.loc[
        (cds_negative_end.iloc[:, GTF_TYPE_IX] == "CDS")
        & (cds_negative_end.iloc[:, GTF_STRAND_IX] == "-")
    ]
    cds_negative_end = (
        cds_negative_end.groupby("transcript_id", group_keys=True)
        .max()
        .iloc[:, GTF_END_IX]
        .reset_index()
        .sort_values("transcript_id")
    )
    cds_negative_end.columns = [cds_positive_start.columns[0], "coord"]
    cds_start = pd.concat([cds_positive_start, cds_negative_start], axis=0).sort_values(
        "transcript_id"
    )

    cds_end = pd.concat([cds_positive_end, cds_negative_end], axis=0).sort_values(
        "transcript_id"
    )

    blast_table = pd.read_csv(blastp_path, sep="\t", header=None)
    blast_table.columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]
    blast_missing_set = np.setdiff1d(
        gtf["transcript_id"].values, np.unique(blast_table.qseqid.values)
    )

    if blast_missing_set.shape[0] >= 1:
        blast_table_missing = pd.DataFrame(
            np.concatenate(
                [
                    np.expand_dims(blast_missing_set, 1),
                    np.full(
                        (blast_missing_set.shape[0], blast_table.shape[1] - 1), np.nan
                    ),
                ],
                axis=1,
            ),
            columns=blast_table.columns,
        )
        blast_table = pd.concat([blast_table, blast_table_missing], axis=0)
    blast_table = blast_table.groupby("qseqid").head(1).sort_values("qseqid")
    protein_fasta = pd.DataFrame(
        {
            "transcript_id": [
                seq_record.id for seq_record in SeqIO.parse(protein_fasta_path, "fasta")
            ],
            "seq": [
                seq_record.seq
                for seq_record in SeqIO.parse(protein_fasta_path, "fasta")
            ],
        }
    ).sort_values("transcript_id")

    pd.DataFrame(
        {
            "Chromosome": gtf.iloc[:, 0].values,
            "Start": gtf.iloc[:, GTF_START_IX].values,
            "Stop": gtf.iloc[:, GTF_END_IX].values,
            "Strand": gtf.iloc[:, GTF_STRAND_IX].values,
            "Source": gtf.iloc[:, GTF_SOURCE_IX].values,
            "CDS_Source": cds_source["source"].values,
            "CDS_Start": cds_start.coord.values,
            "CDS_Stop": cds_end.coord.values,
            "gid": gtf["gene_id"].values,
            "gene_name": gtf["gene_name"].values,
            "tid": orfs["transcript_id"].values,
            "pid": [i.rsplit("|")[0] for i in blast_table["sseqid"].astype(str).values],
            "blastp_identity": blast_table["pident"].values,
            "blastp_bitscore": blast_table["bitscore"].values,
            "transcript_length_nt": orfs["len"].values,
            "orf_length_nt": orfs["orf_len"].values,
            "protein_length_cd": (orfs["orf_len"].values - STOP_CODON_NT)
            / CODON_LENGTH,
            "protein_is_nmd": sqanti_protein.is_nmd.fillna(False).values,
            "protein_splice_category": sqanti_protein.pr_splice_cat.values,
            "protein_splice_subcategory": sqanti_protein.pr_splice_subcat.values,
            "protein_has_stop_codon": orf_completeness.has_stop_codon.values,
            "protein_has_start_codon": orf_completeness.has_start_codon.values,
            "protein_sequence": protein_fasta.seq.values,
        }
    ).to_csv(f"{output_name}", index=False, sep="\t")
    return 0


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--best_orf_path")
    parser.add_argument("--sqanti_protein_path")
    parser.add_argument("--orf_completeness_path")
    parser.add_argument("--output_name")
    parser.add_argument("--gtf_original_path")
    parser.add_argument("--gtf_predicted_path")
    parser.add_argument("--protein_fasta_path")
    parser.add_argument("--blastp_path")

    args = parser.parse_args()
    main(
        best_orf_path=args.best_orf_path,
        sqanti_protein_path=args.sqanti_protein_path,
        orf_completeness_path=args.orf_completeness_path,
        output_name=args.output_name,
        gtf_original_path=args.gtf_original_path,
        gtf_predicted_path=args.gtf_predicted_path,
        protein_fasta_path=args.protein_fasta_path,
        blastp_path=args.blastp_path,
    )
