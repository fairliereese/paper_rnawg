#!/usr/bin/env python
# -*- coding: utf-8 -*-


from argparse import ArgumentParser

import numpy as np
import pandas as pd
from Bio import SeqIO
from typeguard import typechecked

STOP_CODONS = ("TAG", "TAA", "TGA")
START_CODONS = "ATG"


@typechecked
def main(
    input_file_path: str,
    orf_input_seq_path: str,
    output_path: str,
    first_cutoff: float,
    second_cutoff: float,
) -> int:
    full_orf_ids = []
    for record in SeqIO.parse(orf_input_seq_path, "fasta"):
        if str(record.seq).upper().endswith(STOP_CODONS) and str(
            record.seq
        ).upper().startswith(START_CODONS):
            full_orf_ids.append(record.id)
    df = pd.read_csv(f"{input_file_path}", sep="\t")
    df = df.loc[np.isin(df.ID, full_orf_ids)].copy(deep=True)
    df["transcript_id"] = df.ID.str.rsplit("_").str[0]
    df_first_cutoff = df.loc[df["Coding_prob"] > first_cutoff]
    df_second_cutoff = df.loc[df["Coding_prob"] > second_cutoff]
    pd.concat(
        [
            df_first_cutoff.sort_values(["ORF", "Coding_prob"], ascending=False)
            .groupby("transcript_id")
            .head(1),
            df_second_cutoff.sort_values(["ORF", "Coding_prob"], ascending=False)
            .groupby("transcript_id")
            .head(1),
            df.sort_values(["ORF", "Coding_prob"], ascending=False)
            .groupby("transcript_id")
            .head(1),
        ],
        axis=0,
        ignore_index=True,
    ).sort_values(["Coding_prob"], ascending=False).sort_values(
        "Coding_prob", ascending=False
    ).groupby(
        "transcript_id"
    ).head(
        1
    ).to_csv(
        f"{output_path}",
        sep="\t",
        index=False,
        header=True,
    )
    return 0


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input_file_path")
    parser.add_argument("--orf_input_seq_path")
    parser.add_argument("--output_path")
    parser.add_argument("--first_cutoff")
    parser.add_argument("--second_cutoff")

    args = parser.parse_args()
    main(
        input_file_path=args.input_file_path,
        orf_input_seq_path=args.orf_input_seq_path,
        output_path=args.output_path,
        first_cutoff=float(args.first_cutoff),
        second_cutoff=float(args.second_cutoff),
    )
