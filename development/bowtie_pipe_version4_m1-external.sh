#!/bin/bash

# bowtie_pipe_version4_m1 -i e91_and_appris_GRCh38-transc -g ensembl91_and_appris_GRCh38_principal_transc --norc -v 0 -p 4

if [ "$1" = "-i" ]; then
  exp="$2"
  shift 2
else
  { echo "bowtie_pipe needs -i [experiment name] first. not $1. failed" && exit 1; }
fi

if [ "$1" = "-g" ]; then
  index="$2"
  shift 2
else
  { echo "bowtie_pipe needs -g [index name] as second argument. not $1. failed" && exit 1; }
fi

# I am pretty sure I didn't need to do multiple "if" and shift blocks. could've done 1 loop but I don't remember how

for var1
do
   direct=$direct$var1
done

direct=$direct-m1
fprefix="$exp$direct"


index_path="/Volumes/external2/rpf-external/bowtie1-indexes/""$index"
reads="/Volumes/external2/rpf-external/reads/""$exp.fastq"
sam_output="$fprefix.sam" # probably unnecessary variable

mkdir "$index"
cd "$index"

mkdir "$exp-alignment$direct"
cd "$exp-alignment$direct"

echo "experiment: $reads"
echo "alignment index: $index_path"
echo "file prefix: $fprefix"
echo "bowtie alignment command:"
echo "bowtie -S -t $@ -m 1 --max $fprefix-suppressed.fastq $index_path $reads $fprefix.sam"

# redirect stdout/stderr to a file
exec &> logfile.txt

echo "experiment: $reads"
echo "alignment index: $index_path"
echo "file prefix: $fprefix"
echo "bowtie alignment command:"
echo "bowtie -S -t $@ -m 1 --max $fprefix-suppressed.fastq $index_path $reads $fprefix.sam"

bowtie -S -t $@ -m 1 --max "$fprefix-suppressed.fastq" $index_path $reads "$fprefix.sam"

echo "converting to BAM, sorting, indexing"
samtools view -u "$fprefix.sam" | samtools "sort" -o "$fprefix.sorted.bam"
samtools index "$fprefix.sorted.bam"

echo "deleting sam file"
rm "$fprefix.sam"

bed_bam2cov_tsv.sh -bam "$fprefix.sorted.bam" -bed "hsp90s.bed" -cov "$fprefix.tsv.gz"
cov_tsvgz2csv.py "$fprefix.tsv.gz"
samtools idxstats "$fprefix.sorted.bam" > "$fprefix-aligns_per_gene.txt"

echo "aligning reads suppressed my m=1 (mulimappers)"
echo "reporting all alignments"

echo "bowtie -S -t $@ -a $index_path ./$fprefix-suppressed.fastq $fprefix-a_suppr_realigned.sam"
bowtie -S -t $@ -a $index_path "./$fprefix-suppressed.fastq" "$fprefix-a_suppr_realigned.sam"

echo "converting to BAM, sorting, indexing"
samtools view -u "$fprefix-a_suppr_realigned.sam" | samtools "sort" -o "$fprefix-a_suppr_realigned.sorted.bam"
samtools index "$fprefix-a_suppr_realigned.sorted.bam"

echo "deleting sam file"
rm "$fprefix-a_suppr_realigned.sam"

samtools idxstats "$fprefix-a_suppr_realigned.sorted.bam" > "$fprefix-a_suppr_realigned-aligns_per_gene.txt"