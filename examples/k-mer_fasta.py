#!/usr/bin/env python
# change the shebang to your specific python3 binary location if this doesn't work

'''
Usage:
k-mer_fasta.py [k] [input name] [output name]
'''

from Bio import SeqIO
import sys


def gen_reads(seq, k):
    '''
    gen_reads generates list of k base "reads" sliding across transcripts.
    '''
    k2 = k - 1
    reads = [''] * (len(seq) - k2)
    for i in range(len(seq) - k2):
        reads[i] = seq[i:i + k]
    return reads


def read_to_fasta(reads, transcript_name, k):
    '''
    read_to_fasta converts reads to fasta format for each transcript
    '''
    fasta = [''] * (len(reads) * 2)
    j = 0
    for n, s in enumerate(reads):
        # info listed before each read.
        # position of alignment in SAM file is 1-based. So position is n+1
        fasta[j] = '>{}.{}'.format(transcript_name, str(n + 1))
        fasta[j + 1] = s
        j = j + 2
    return fasta


def main():
    k = int(sys.argv[1])
    # import and parse fasta file with transcripts. returns seqs: list of sequences
    seqs = []
    ids = []
    transcriptome_file = str(sys.argv[2])

    with open(transcriptome_file) as handle:
        fasta_sequences = SeqIO.parse(handle, 'fasta')
        for fasta in fasta_sequences:
            seqs.append(str(fasta.seq))
            ids.append(str(fasta.id))

    with open(str(sys.argv[3]), 'w+') as output_handle:
        for sequence, name in zip(seqs, ids):
            fasta_reads = read_to_fasta(gen_reads(sequence, k), name, k)
            for i in fasta_reads:
                output_handle.write("{}\n".format(i))


if __name__ == "__main__":
    main()
