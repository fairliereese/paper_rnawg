#!/usr/bin/env python
# -*- coding: utf-8 -*-


from argparse import ArgumentParser

import numpy as np
import pandas as pd
from Bio import SeqIO
from typeguard import typechecked

GTF_GENE_INFO_IX = 8

STOP_CODONS = ("TAG", "TAA", "TGA")
START_CODONS = "ATG"


@typechecked
def main(
    cpat_seqs: str,
    orfanage_seqs: str,
    cpat_info: str,
    orfanage_info: str,
    output_path: str,
) -> int:
    df_cpat = pd.read_csv(cpat_info, sep="\t")
    df_orfanage = pd.read_csv(orfanage_info, sep="\t", header=None)
    df_orfanage["transcript_id"] = (
        df_orfanage.iloc[:, GTF_GENE_INFO_IX]
        .str.rsplit(";")
        .str[0]
        .str.rsplit('"')
        .str[1]
        .values
    )
    orfanage_cds = np.unique(df_orfanage.transcript_id)
    transcript_id = []
    source = []
    has_start_codon = []
    has_stop_codon = []

    for record in SeqIO.parse(orfanage_seqs, "fasta"):
        if record.id in orfanage_cds:
            transcript_id.append(record.id)
            source.append("ORFanage")
            has_stop_codon.append(
                str(record.seq).upper().endswith(STOP_CODONS)
            )
            has_start_codon.append(
                str(record.seq).upper().startswith(START_CODONS)
            )
    cpat_cds = np.unique(df_cpat.ID.values)
    for record in SeqIO.parse(cpat_seqs, "fasta"):
        if record.id in cpat_cds:
            transcript_id.append(record.id.rsplit("_")[0])
            source.append("CPAT")
            has_stop_codon.append(
                str(record.seq).upper().endswith(STOP_CODONS)
            )
            has_start_codon.append(
                str(record.seq).upper().startswith(START_CODONS)
            )

    pd.DataFrame(
        {
            "transcript_id": transcript_id,
            "source": source,
            "has_stop_codon": has_stop_codon,
            "has_start_codon": has_start_codon,
        }
    ).to_csv(f"{output_path}", sep="\t", index=False, header=True)
    return 0


if __name__ == "__main__":

    parser = ArgumentParser()

    parser.add_argument("--cpat_seqs")
    parser.add_argument("--orfanage_seqs")
    parser.add_argument("--cpat_info")
    parser.add_argument("--orfanage_info")
    parser.add_argument("--output_path")

    args = parser.parse_args()
    main(
        cpat_seqs=args.cpat_seqs,
        orfanage_seqs=args.orfanage_seqs,
        cpat_info=args.cpat_info,
        orfanage_info=args.orfanage_info,
        output_path=args.output_path,
    )
