'''
Take known list of multi-mapping transcripts from yeast transcriptome (previous analysis)
identify the multi-mapping transcripts in the example transcriptome
'''


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob, os

# %%
# ==============================================================================
# ==============================================================================

fasta_path = './example_transcriptome.fasta'
ids = []
from Bio import SeqIO
with open(fasta_path) as handle:
    fasta_sequences = SeqIO.parse(handle, 'fasta')
    for fasta in fasta_sequences:
        ids.append(str(fasta.id))

mms = pd.read_csv('MMtable_merged_30mers_v2_version1.csv')
mms = list(mms['transcript'])



# %%
# ==============================================================================
# ==============================================================================
mmers = []
for i in ids:
    if i in mms:
        mmers.append(i)
len(mmers)
len(mms)
len(ids)


