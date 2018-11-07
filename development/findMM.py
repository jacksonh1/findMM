#!/usr/bin/env python

#./argparse_test2.py -p4 -i ./index/example_transcriptome -k 30 -ref example_transcriptome.fasta

# Parsing command line arguments
import argparse
import sys
import os
import subprocess


parser = argparse.ArgumentParser(
    description='''identify multi-mapping transcripts and their connectivity in a reference transcriptome fasta file'''
)
parser.add_argument(
    '-b',
    metavar="<bowtie path>",
    default='bowtie',
    help='''the path to the bowtie executable. Default: bowtie in PATH'''
)
parser.add_argument(
    '-p',
    '--threads',
    metavar="<INT>",
    default=1,
    type=int,
    help='''Number of alignment threads for bowtie to use (given directly to bowtie parameter: -p). Default: 1'''
)
parser.add_argument(
    '-out',
    metavar='<directory>',
    default='./',
    help='''Directory where script outputs are saved. Default is current directory'''
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
    help='''k-mer read lengths (i.e. k). Recommend k be equal to the smallest read length in an experimental dataset'''
)
parser.add_argument(
    '-ref',
    '--reference',
    metavar='<file>',
    required=True,
    help='''input reference transcriptome file in fasta format'''
)

args = parser.parse_args()

print('''
Running findMM with the following parameters:
 - bowtie path: {}
 - threads: {}
 - output directory: {}
 - bowtie index: {}
 - k: {}
 - reference transcriptome: {}
'''.format(args.b, args.threads, args.out, args.index, args.k, args.reference))


# %%
# ==============================================================================
# // run scripts
# ==============================================================================

def k_mer_generation(args):
    # parse input fasta file name to get output filename
    base = os.path.basename(args.reference)
    basenoext, ext = os.path.splitext(base)
    # could put extension check here
    directory = os.path.dirname(args.reference)
    new_filename = '{}-{}-mers{}'.format(basenoext, str(args.k), '.fa')
    k_mer_file = os.path.join(args.out, new_filename)
    subprocess.call(
        '../scripts/k-mer_fasta.py {} {} {}'.format(
            args.k,
            args.reference,
            k_mer_file
        ),
        shell=True
    )
    return k_mer_file


print('generating transcriptome k-mers')
k_mer_file = k_mer_generation(args)

# %%
# ==============================================================================
# // TITLE
# ==============================================================================
# eggs = 'eggs2'
#
#
# try:
#     if egg != 'eggs':
#         raise ValueError(
#             'input reference file does not appear to be a "fasta" file')
# except ValueError as er:
#     print(er)
#     check = input('Would you like to continue anyway? (Y or N)'
