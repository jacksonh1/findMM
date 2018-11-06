'''
- k
- mismatches
    - default: `-v 2`
- bowtie path
    - default: uses bowtie in PATH
- threads
    - default: 2?
- index
	- default: make one
- input transcriptome fasta
- output_path
    - k-mer reads go here
    - table and network go here
    - default: cwd
'''
import argparse
import sys

parser = argparse.ArgumentParser(description='identify multi-mapping transcripts and their connectivity in a reference transcriptome fasta file')
parser.add_argument('k', metavar='k', nargs='+', type=int, help='k-mer read lengths (i.e. k). Recommend k be equal to the smallest read length in an experimental dataset')
parser.add_argument('bowtie path', default='bowtie', type=str, help='the path to the bowtie executable. Default is bowtie in PATH')
args = parser.parse_args()