# findMM

run `./findMM.py -h` for info and command line arguments <br>
See `example` folder for an example use of findMM <br>

## overview
findMM is be a simple, python-based command line utility for identifying the multi-mapping transcripts and their multi-mapping connectivity to other transcripts in a reference transcriptome fasta file. The approach to identifying the multi-mapping transcripts is similar to the [`crossmap`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.crossmap.html) script from the plastid package.

findMM takes a reference transcriptome as input and produces a table of all of the multi-mapping transcripts. For each transcript the multi-mapping table includes the transcript name, transcript length, how many other transcripts it multi-maps with, if there are any internal multi-maps, and what percent of the transcript sequence multi-maps.
<br>Multi-mapping table:

| transcript | length (bp) | percent of transcript that multi-maps | external multi-maps | internal multi-maps? |
| ---------- | ----------- | ------------------------------------ | ------------------ | -------------------- |
| YJL138C  	 | 1188        | 100                                  | 3                  | no                   |
| ...      	 | ...         | ...                                  | ...                | ...                  |

findMM also outputs a multi-mapping network file describing the multi-mapping patterns within the transcriptome. The network file is simply a `csv` file where each row is a multi-mapping connection including both transcript names and the number of k-mers which map to both transcripts. The network can be visualized in [cytoscape](https://cytoscape.org/).
<br>Multi-mapping network:

| read_origin | transcript | MM-kmer-alignments |
| ----------  | ---------- | ------------------ |
| YLL024C  	  | YAL005C    | 1567               |
| ...      	  | ...        | ...                |
- **read_origin** - transcript the k-mers originate from
- **transcript** - transcript the k-mers multi-map with
- **MM-kmer-alignments** - Number of k-mers that multi-map between the two transcripts


## Prerequisites for findMM:
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
- [Samtools](http://www.htslib.org/)
- python version 3.6 or later
- the following python modules:
	- [biopython](https://biopython.org/wiki/Download)
	- [pandas](https://pandas.pydata.org/pandas-docs/stable/install.html)
- samtools and python available in PATH variable

## using findMM:

To run fingMM, use the script `findMM.py`. You can add `findMM.py` to your PATH for convenience.

For example, with `findMM.py` in current directory, you would run:
```bash
./findMM.py [options]
```
where `[options]` are the command line arguments<br>
Running `./findMM.py -h` lists the available command line arguments and explains what each parameter is:
```bash
./findMM.py -h
```
```
usage: findMM.py [-h] -ref <file> -k <INT> -i <ebwt_base> [-b <bowtie path>]
                 [-v <INT>] [-p <INT>] [-out <directory>]

identify multi-mapping transcripts and their connectivity in a reference
        transcriptome fasta file

optional arguments:
  -h, --help            show this help message and exit
  -ref <file>, --reference <file>
                        input reference transcriptome file in fasta format
  -k <INT>              k-mer read lengths (i.e. k). Recommend k be equal to
                        the smallest read length in an experimental dataset
  -i <ebwt_base>, --index <ebwt_base>
                        location of bowtie index made from input reference transcriptome.
                        To create bowtie index you can use:
                        bowtie-build <reference_in> <ebwt_outfile_base>
  -b <bowtie path>      the path to the bowtie executable. Default: bowtie in
                        PATH
  -v <INT>, --mismatches <INT>
                        maximum number of mismatches allowed in alignments
                        (given directly to bowtie parameter -v). Default: 2
  -p <INT>, --threads <INT>
                        Number of alignment threads for bowtie to use
                        (given directly to bowtie parameter -p). Default: 1
  -out <directory>      Directory where script outputs are saved. Default is
                        current directory
```
**The following arguments are required:**
- `-ref <file>, --reference <file>`
- `-k <INT>`
- `-i <ebwt_base>, --index <ebwt_base>`


