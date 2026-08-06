# copy-number-prediction
Python script to predict gene copy number in prokaryotic genome.

Set up a conda environment with the following dependencies installed.

## Dependencies:
- bedtools version >=2.31.1
- bwa version >=0.7.19
- samtools >=1.23

## Usage
python combined_run.py \
    -f sample.fa \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -g sample.gff \
    -o sample.tsv \
    -q <query_gene_name> \
    -s <subject_gene_name>

## Parameters
```
  ### -f <sample.fa> your genome assembly FASTA
  ### -1 <sample_R1.fastq.gz> read R1
  ### -2 <sample_R2.fastq.gz> read R2
  ### -g <sample.gff> annotation file from Prokka
  ### -o <sample.tsv> name of output file
  ### -q <gene name> Query gene, usually single-copy gene to use as baseline
  ### -s <gene name> Subject gene, to calculate copy number
```
