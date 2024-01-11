#!/usr/bin/env python3

# make GTF file that includes the ORF regions (as CDS features)
# input - pacbio gtf ('jurkat.collapsed.gtf'), orf calls ('jurkat_refine_orf_calls.tsv')
# output - pacbio gtf with added "cds" features (orfs)

import argparse
import copy
import warnings
from collections import defaultdict

import gtfparse
import numpy as np
import pandas as pd


def string_to_boolean(string):
    """
    Converts string to boolean

    Parameters
    ----------
    string :str
    input string

    Returns
    ----------
    result : bool
    output boolean
    """
    if isinstance(string, bool):
        return str
    if string.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif string.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def get_first_block_index(orf_coord, cblens, pblens):
    # get the index corresponding to the first block containing the orf start
    # return index, and the dela (spacing upstream of end)
    for i, cblen in enumerate(cblens):
        if orf_coord <= cblen:
            delta = cblen - orf_coord
            return i, delta
    return i, 0


def make_cds_coords_positive_strand(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1 : i2 + 1]
    # trim ends to orf start/end
    orf_coords[0][0] = orf_coords[0][1] - delta1
    orf_coords[-1][1] = orf_coords[-1][1] - delta2
    return orf_coords


def make_cds_coords_negative_strand(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1 : i2 + 1]
    # trim ends to orf start/end
    orf_coords[0][1] = orf_coords[0][0] + delta1
    orf_coords[-1][0] = orf_coords[-1][0] + delta2
    return orf_coords


def get_min_and_max_coords_from_exon_chain(coords):
    """Gets the minumum and maximum coordinates from exon chain

    Args:
        coords ([[int, int]]): [[start,end]] exon coordinate chain

    Returns:
        (int, int): start and end coordinates of entire chain
    """

    min_coord = min(map(min, coords))
    max_coord = max(map(max, coords))
    return min_coord, max_coord


def make_pacbio_cds_gtf(sample_gtf, called_orfs, name):
    """Makes PacBio CDS and saves file with CDS

    Args:
        sample_gtf (filename): sample_gtf file
        refined_orfs (filename): aggregate_orf info. from Refined DB
        called_orfs (filename): orf calls from ORF_Calling
        pb_gene (filename): PacBio gene name cross reference
        name (string): name of sample
    """
    # import gtf, only exon info.
    # only move forward with representative pb isoform (for same-protein groups)
    gtf = gtfparse.read_gtf(sample_gtf)
    gtf_gene_mapping = gtf.loc[gtf["feature"] == "transcript"]
    gtf_gene_mapping = gtf_gene_mapping[["transcript_id", "gene_id"]]
    gtf = gtf[
        ["seqname", "feature", "start", "end", "strand", "transcript_id"]
    ]
    gtf = gtf.loc[gtf["feature"] == "exon"]
    gtf.columns = ["chr", "feat", "start", "end", "strand", "acc"]
    # only move forward with "base accession" (representative pb)
    # pb coords into dict
    pbs = defaultdict(
        lambda: ["chr", "strand", [], [], [], []]
    )  # pb -> [chr, strand, [start, end], [block lengths],[cum. block lengths], [prior cumulative block lengths]]
    # PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
    for i, row in gtf.iterrows():
        chr, feat, start, end, strand, acc = row
        pbs[acc][0] = chr
        pbs[acc][1] = strand
        pbs[acc][2].append([int(start), int(end)])
    # sort all coords, calc blocks
    for acc, infos in pbs.items():
        strand = infos[1]
        if strand == "+":
            infos[2] = sorted(infos[2])
        elif strand == "-":
            infos[2] = sorted(infos[2], reverse=True)
        infos[3] = np.array([end - start + 1 for [start, end] in infos[2]])
        infos[4] = np.cumsum(infos[3])
        infos[5] = infos[4] - infos[3]

    # read in the ranges of orf on pb transcripts
    ranges = pd.read_table(called_orfs)[
        ["transcript_id", "ORF_start", "ORF_end"]
    ]
    # print(gtf_gene_mapping)
    # print(gtf_gene_mapping.gene_id)
    # pb_gene = pd.Series(
    #    gtf_gene_mapping["gene_id"], index=gtf_gene_mapping["transcript_id"]
    # ).to_dict()
    pb_gene = dict(
        zip(gtf_gene_mapping["transcript_id"], gtf_gene_mapping["gene_id"])
    )
    with open(f"{name}_cpat_with_cds.gtf", "w") as ofile:
        for i, row in ranges.iterrows():
            acc, orf_start, orf_end = row
            # remove stop exon
            orf_end = orf_end - 3
            if acc in pbs:
                if acc in pb_gene.keys():
                    gene = pb_gene[acc]
                else:
                    raise ValueError
                    gene = "-"
                infos = pbs[acc]
                chr, strand, coords, blens, cblens, pblens = infos

                i1, delta1 = get_first_block_index(orf_start, cblens, pblens)
                i2, delta2 = get_first_block_index(orf_end, cblens, pblens)
                if strand == "+":
                    orf_coords = make_cds_coords_positive_strand(
                        i1, delta1, i2, delta2, coords
                    )
                elif strand == "-":
                    orf_coords = make_cds_coords_negative_strand(
                        i1, delta1, i2, delta2, coords
                    )
                # write out the coordinates
                # out_acc = f'transcript_id "{acc}";'
                # acc_w_gene_w_cpm = gene + "|" + acc
                out_acc = f'transcript_id "{acc}"; gene_id "{gene}";'

                out_acc_exon = f'transcript_id "{acc}";'

                tstart, tend = get_min_and_max_coords_from_exon_chain(coords)
                ofile.write(
                    "\t".join(
                        [
                            chr,
                            "hg38_canon",
                            "transcript",
                            str(tstart),
                            str(tend),
                            ".",
                            strand,
                            ".",
                            out_acc,
                        ]
                    )
                    + "\n"
                )
                for [start, end] in coords:
                    ofile.write(
                        "\t".join(
                            [
                                chr,
                                "hg38_canon",
                                "exon",
                                str(start),
                                str(end),
                                ".",
                                strand,
                                ".",
                                out_acc_exon,
                            ]
                        )
                        + "\n"
                    )
                for [start, end] in orf_coords:
                    ofile.write(
                        "\t".join(
                            [
                                chr,
                                "CPAT",
                                "CDS",
                                str(start),
                                str(end),
                                ".",
                                strand,
                                ".",
                                out_acc_exon,
                            ]
                        )
                        + "\n"
                    )


def main():
    parser = argparse.ArgumentParser(
        "IO file locations for make pacbio cds gtf"
    )
    parser.add_argument(
        "--name",
        action="store",
        dest="name",
        help="name of sample - used for output file name",
    )
    parser.add_argument(
        "--sample_gtf",
        action="store",
        dest="sample_gtf",
        help="sample gtf, from sqanti3",
    )
    parser.add_argument(
        "--called_orfs",
        action="store",
        dest="called_orfs",
        help="agg orf tsv file, from refined DB",
    )

    results = parser.parse_args()
    warnings.simplefilter(action="ignore", category=FutureWarning)
    make_pacbio_cds_gtf(
        results.sample_gtf,
        results.called_orfs,
        results.name,
    )


if __name__ == "__main__":
    main()
