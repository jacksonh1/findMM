# findMM

## this project is not yet ready for distribution, however the script: `findMM.py` works
run `./findMM.py -h` for info and command line arguments <br>
See `example` folder for an example use of findMM <br>
**all suggestions and constructive feedback are welcome and appreciated**

## overview
---
findMM is be a simple, python-based command line utility for identifying the multi-mapping transcripts and their multi-mapping connectivity to other transcripts in a reference transcriptome fasta file. The approach to identifying the multi-mapping transcripts is similar to the [`crossmap`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.crossmap.html) script from the plastid package.

findMM takes a reference transcriptome as input and outputs a table of all of the multi-mapping transcripts in addition to a multi-mapping network `csv` file describing the multi-mapping patterns within the transcriptome. For each transcript the multi-mapping table includes how many other transcripts it multi-maps with, if there are any internal multi-maps, and what percent of the transcript sequence multi-maps.

## Prerequisites for findMM:
---
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
- [Samtools](http://www.htslib.org/)
- python version 3.6 or later
- the following python modules:
	- [biopython](https://biopython.org/wiki/Download)
	- [pandas](https://pandas.pydata.org/pandas-docs/stable/install.html)
- samtools and python available in PATH variable
---



### major/essential TO-DOs
- clean up generated files if script fails


