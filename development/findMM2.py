#!/usr/bin/env python

#./findMM.py -p4 -i ./index/example_transcriptome -k 30 -ref example_transcriptome.fasta

# could add 'save_kmers'. Action.... In kmer_function(args, save_kmers):

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
    '-v',
    '--mismatches',
    metavar="<INT>",
    default=2,
    type=int,
    help='''maximum number of mismatches allowed in alignments (given directly to bowtie parameter -v). Default: 2'''
)
parser.add_argument(
    '-p',
    '--threads',
    metavar="<INT>",
    default=1,
    type=int,
    help='''Number of alignment threads for bowtie to use (given directly to bowtie parameter -p). Default: 1'''
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
 - number of allowed mismatches (-v): {}
 - threads used in alignment: {}
 - output directory: {}
 - bowtie index: {}
 - k: {}
 - reference transcriptome: {}
'''.format(args.b, args.mismatches, args.threads, args.out, args.index, args.k, args.reference))

dict(args._get_kwargs())
# %%
# ==============================================================================
# // run scripts
# ==============================================================================


def k_mer_generation(args):
    subprocess.call('mkdir {}'.format(args.out), shell=True)
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


def bowtie_alignment1(args, k_mer_file):
    '''
    args - argparse parsed command line arguments
    k_mer_file - path to fasta formatted k-mer reads
    '''
    from string import Template
    params = dict(args._get_kwargs())
    # subprocess.call('mkdir {}'.format(args.out), shell=True)
    k_mer_noext, ext = os.path.splitext(k_mer_file)
    sam_path = '{}.sam'.format(k_mer_noext)
    max_output = '{}-multi-mapping_kmers.fastq'.format(k_mer_noext)
    sorted_bam_path = '{}_sorted.bam'.format(k_mer_noext)
    MM_sam_path = '{}_MMers.sam'.format(
        k_mer_noext)  # MMing k_mer alignment sam
    params['sam_path'] = sam_path
    params['max_output'] = max_output
    params['k_mer_file'] = k_mer_file
    params['sorted_bam_path'] = sorted_bam_path
    params['MM_sam_path'] = MM_sam_path
    bowtie_command1 = Template(
        '$b -S -p $threads --norc -v $mismatches -f -m 1 --max $max_output $index $k_mer_file $sam_path')
    print('running bowtie alignment:\n' + bowtie_command1.substitute(params))
    subprocess.call(
        bowtie_command1.substitute(params),
        shell=True
    )
    # TODO - change to log file and display
    # print('converting to BAM, sorting, and indexing')
    # subprocess.call(
    #     'samtools view -u {} | samtools sort -o {}'.format(
    #         params['sam_path'], params['sorted_bam_path']),
    #     shell=True
    # )
    # subprocess.call('samtools index {}'.format(
    #     params['sorted_bam_path']), shell=True)
    subprocess.call('rm {}'.format(params['sam_path']), shell=True)

    # bowtie alignment 2
    # MM_sorted_sam_path = '{}_MMers_sorted.sam'.format(k_mer_noext)# MMing k_mer alignment BAM
    # params['MM_sorted_sam_path'] = MM_sorted_sam_path
    bowtie_command2 = Template(
        '$b -S -p $threads --norc -v $mismatches -f -a $index $max_output $MM_sam_path')  # realign just multi-mapping k-mers
    print('running bowtie alignment with multi-mapping k-mers: \n' +
          bowtie_command2.substitute(params))
    subprocess.call(
        bowtie_command2.substitute(params),
        shell=True
    )
    return params


def main():
    print('generating transcriptome k-mers')
    k_mer_file = k_mer_generation(args)
    params = bowtie_alignment1(args, k_mer_file)


if __name__ == "__main__":
    main()
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
