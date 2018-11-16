#!/bin/bash

usage: findMM.py [-h] [-b <bowtie path>] [-p <INT>] [-out <directory>] -i
                 <ebwt_base> -k <INT> -ref <file>

identify multi-mapping transcripts and their connectivity in a reference
transcriptome fasta file

optional arguments:
  -h, --help            show this help message and exit
  -b <bowtie path>      the path to the bowtie executable. Default: bowtie in
                        PATH
  -p <INT>, --threads <INT>
                        Number of alignment threads for bowtie to use (given
                        directly to bowtie parameter: -p). Default: 1
  -out <directory>      Directory where script outputs are saved. Default is
                        current directory
  -i <ebwt_base>, --index <ebwt_base>
                        location of bowtie index made from input reference
                        transcriptome. Use bowtie-build <reference_in>
                        <ebwt_outfile_base> example: bowtie-build
                        ./example_transcriptome.fasta
                        ./index/example_transcriptome future: make argument
                        optional and build index in ./index by default
  -k <INT>              k-mer read lengths (i.e. k). Recommend k be equal to
                        the smallest read length in an experimental dataset
  -ref <file>, --reference <file>
                        input reference transcriptome file in fasta format

```bash
mkdir [output directory]
bowtie -S -p [] -f -m 1 --max [multi-mapping_kmers.fastq] [index_path] [k_mer_file] [output directory + output_sam_file]
```

inputs:
`threads`
`--max output path = [k_mer_name + multi-mapping_kmers.fastq]`
`index_path`
`k_mer_file`
sam path = `[output directory + output_sam_file]`
