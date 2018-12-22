# findMM

findMM is be a simple, python-based command line utility for identifying the multi-mapping transcripts and their multi-mapping connectivity to other transcripts in a reference transcriptome fasta file. The approach to identifying the multi-mapping transcripts is similar to the [`crossmap`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.crossmap.html) script from the plastid package.

findMM takes a reference transcriptome as input and outputs a table of all of the multi-mapping transcripts in addition to a multi-mapping network `csv` file describing the multi-mapping patterns within the transcriptome. For each transcript the multi-mapping table includes how many other transcripts it multi-maps with, if there are any internal multi-maps, and what percent of the transcript sequence multi-maps.

Requirements:
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
- [Samtools](http://www.htslib.org/)
- python 3.6 or later
	- [biopython](https://biopython.org/wiki/Download)
	- [pandas](https://pandas.pydata.org/pandas-docs/stable/install.html)
- bowtie, samtools, and python in PATH

---

## this project is not yet ready for distribution, however the script: findMM.py works
**all suggestions and constructive feedback are welcome and appreciated**

### major/essential TO-DOs
- clean up generated files if script fails

	
