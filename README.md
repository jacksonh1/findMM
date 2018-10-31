# findMM

findMM will be a simple, python-based command line utility for identifying the multi-mapping transcripts and their multi-mapping connectivity to other transcripts in a reference transcriptome fasta file. The approach to identifying the multi-mapping transcripts is similar to the [`crossmap`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.crossmap.html) script from the plastid package. I highly recommend `crossmap` for working with a reference genome. Also if you want a mask file specifying the specific multi-mapping regions of the genome.

findMM takes a reference transcriptome as input and outputs a table of all of the multi-mapping transcripts. For each transcripts it includes how many other transcripts it multi-maps with and if there are any internal multi-maps.

Example output:
multi-mapping transcript | number of external multi-maps | internal multi-maps
--- | --- | ---
example | 5 | yes

Requirements:
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
- [Samtools](http://www.htslib.org/)
- python 3.6 or later
- bowtie, samtools, and python in PATH
