# copy-number-prediction
Python script to predict gene copy number in prokaryotic genome.

Set up a conda environment with the following dependencies installed.

## Dependencies:
- bedtools >=2.31.1
- bwa >=0.7.19
- samtools >=1.23

## Input
The tool takes the following as input:
- your genome FASTA file
- your Illumina PE reads
- your genome annotation file (.gff/.gff3)


## Usage
```
% python copy_number.py \
    -f sample.fa \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -g sample.gff \
    -o sample.tsv \
    -q <query_gene_name> \
    -s <subject_gene_name>
    -t <number_of_threads>

OPTIONS:
  -f <sample.fa> your genome assembly from Prokka/Bakta
  -1 <sample_R1.fastq.gz> read R1
  -2 <sample_R2.fastq.gz> read R2
  -g <sample.gff> annotation file from Prokka/Bakta
  -o <sample.tsv> name of output file
  -q <gene name> Query gene, usually single-copy gene to use as baseline
  -s <gene name> Subject gene, to calculate copy number
  -t <number_of_threads> Set number of threads to run
```

## Output
The tool produces a summary output TSV which contains information on the depth of your reference gene, depth of your gene of interest, and the predicted copy number of your gene of interest.

```
sample	query_gene	query_depth	subject_contig	subject_gene	subject_depth	copy_number
AUSMDU00075849	rpoB	399.33	contig_6	gyrA	292.93	0.734
```
