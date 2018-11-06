# code refinement and testing for distribution
---

`example_transcriptome.fasta` is just a smaller transcriptome created for development purposes


**plan brainstorming:**
- use `argparse` for CLI
- generate k-mers and output k-mer read fasta file
- use bowtie for multi-mapper identification
- build multi-mapping network for generating table
    - parse SAM file
    - probably have to use pandas
- build table
    - parse network
- output table - csv
- output multi-mapping network

**inputs**:
- k
- output_path
    - k-mer reads go here
    - table and network go here
    - default: cwd
- mismatches
    - default: `-v 2`
- bowtie path
    - default: uses bowtie in PATH
- threads
    - default: 2?
- index
    - default: bowtie-build?
        - `./index/`

Documentation:
- python path
    - default uses `python` in PATH
    - change shebang if you use a different python





bowtie alignment test: `bowtie -S -t --norc -v 2 -p 4 -f -m 1 ./index/example_transcriptome ./30-mers_example.fa 30-mers_example.sam`
