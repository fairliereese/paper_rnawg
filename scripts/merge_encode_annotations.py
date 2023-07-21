#!/usr/bin/python3
from argparse import ArgumentParser
import os
import sys
from xopen import xopen


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.output:
        target = xopen(args.output, 'wt')
    else:
        target = sys.stdout

    for filename in args.filenames:
        filetype = detect_file_type(filename)
        with xopen(filename, 'rt') as instream:
            if filetype == 'unrecognized':
                parser.error("Unrecognized file type")
            elif filetype == 'gtf':
                cat(instream, target)
            elif filetype == 'fasta':
                fasta_cat(instream, target)

    if args.output:
        target.close()


def make_parser():
    parser = ArgumentParser(
        description="Utility to combine multiple annotation GTFs and "
                    "synthesize GTF records for spike-ins from fasta files"
    )
    parser.add_argument('filenames', nargs='+', help='name of GTF and fasta files to combine.')
    parser.add_argument('-o', '--output', help='Destination for combined file')
    return parser


def detect_file_type(filename):
    """Read a file and guess if its a fasta or gtf file
    See detect_file_type_stream for details
    """
    stream = xopen(filename, 'rt')
    filetype = detect_file_type_stream(stream)
    stream.close()
    return filetype


def detect_file_type_stream(stream):
    """Read a file and guess if its a fasta or gtf file
    Parameters
    ----------
    stream: file
        name of file object to examine
    Returns
    -------
    str
        fasta, gtf, or unrecognized
    """
    head = stream.read(1024)
    if head[0] == '>':
        return 'fasta'
    elif '\t' in head:
        return 'gtf'
    else:
        return 'unrecognized'


def cat(instream, outstream, strip_comments=True):
    """Copy a file from one stream to another
    """
    for line in instream:
        if not (strip_comments and line.startswith('#')):
            outstream.write(line)


def fasta_cat(instream, outstream):
    """Synthesize GTF records from a fasta stream
    """
    name = None
    seq_list = []
    for line in instream:
        line = line.strip()
        if line[0] == '>':
            if len(seq_list) > 0:
                write_fasta_gtf(outstream, name, seq_list)
                seq_list = []
            name = line[1:]
        else:
            seq_list.append(line)

    write_fasta_gtf(outstream, name, seq_list)


def write_fasta_gtf(outstream, name, seq_list):
    seq = ''.join(seq_list)

    attributes = ' '.join([
        'gene_id "gSpikein_{}";'.format(name),
        'transcript_id "tSpikein_{}";'.format(name),
    ])
    outstream.write('\t'.join([
        name,
        'spikein',
        'exon',
        '1',
        str(len(seq)),
        '.',
        '+',
        '.',
        attributes
    ]))
    outstream.write(os.linesep)


if __name__ == '__main__':
    main()
