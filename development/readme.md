# code refinement and testing for distribution
---

`example_transcriptome.fasta` is just a smaller transcriptome created for development purposes

Question: how do I deal with permissions? Each script probably won't be executable by default.
Check out other programs to see how they do it

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
- k - **required**
- mismatches - **optional**
    - default: `-v 2`
- bowtie path - **optional**
    - default: uses bowtie in PATH
- threads - **optional**
    - default: 2?
- index - **optional**
    - default: bowtie-build?
        - `./index/`
- input transcriptome fasta file - **required**
    - positional
- output_path - **optional**
    - k-mer reads go here
    - table and network go here
    - default: cwd

**Documentation**:
- python path
    - default uses `python` in PATH
    - change shebang if you use a different python


major steps in script:
1. generate k_mers
    - inputs
        - reference transcriptome
        - k
    - output
        - k-mer fasta file
2. optional - generate index
    - inputs
        - reference transcriptome
    - output
        - bowtie index
3. align back to transcriptome with `-m 1 --max`
    - input
        - k-mer fasta dataset
        - reference transcriptome index
    - output
        - multi-mapping k-mers as fasta file
            - could generate simple list from this
        - alignment sam file
4. align multi-mappers back to transcriptome with `-a`
    - input
        - multi-mapping k-mers suppressed by `m=1` in last step
        - reference transcriptome index
    - output
        - alignment sam file
5. network construction
    - input
        - alignment sam file from 4
    - output
        - more detailed table `csv`
        - network `csv`




bowtie alignment test: `bowtie -S -t --norc -v 2 -p 4 -f -m 1 ./index/example_transcriptome ./30-mers_example.fa 30-mers_example.sam`



## first bowtie alignment:

inputs:
threads
max_output = [k_mer_name + multi-mapping_kmers.fastq]`
index_path
k_mer_file
sam_path = `[output directory + output_sam_file]`
norc
mismatches

```bash
mkdir [output directory]
$bowtie -S -p $threads --norc $mismatches -f -m 1 --max $max_output $index_path $k_mer_file $sam_path
```
