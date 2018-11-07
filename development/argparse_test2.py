#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(
    description='''identify multi-mapping transcripts and their connectivity in a reference transcriptome fasta file'''
)
parser.add_argument(
    '-b',
    metavar="<bowtie path>",
    default='bowtie',
    help='the path to the bowtie executable. Default: bowtie in PATH'
)
parser.add_argument(
    '-p',
    '--threads',
    default=1,
    type=int,
    help='Number of alignment threads for bowtie to use (given directly to bowtie parameter: -p). Default: 1'
)
parser.add_argument(
    '-out',
    metavar='<directory>',
    default='./',
    help='directory where script outputs are saved. Default is current directory'
)
parser.add_argument(
    '-i',
    '--index',
    required=True,
    metavar="<ebwt_base>",
    help='''location of bowtie index made from input reference transcriptome.
    Use bowtie-build <reference_in> <ebwt_outfile_base>
    example:
    bowtie-build ./example_transcriptome.fasta ./index/example_transcriptome
    future:
    make argument optional and build index in ./index by default'''
)
parser.add_argument(
    '-k',
    metavar='<INT>',
    required=True,
    type=int,
    help='k-mer read lengths (i.e. k). Recommend k be equal to the smallest read length in an experimental dataset'
)
parser.add_argument(
    '-in',
    metavar='<file>',
    required=True,
    help='input reference transcriptome file in fasta format'
)

args = parser.parse_args()

# print(args.k, args)
