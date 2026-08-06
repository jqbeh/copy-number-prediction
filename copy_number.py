#!/usr/bin/env python3

"""
Usage:

python copy_number.py \
    -f sample.fa \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -g sample.gff \
    -o output.tsv \
    -q rpoB \
    -s OXA-23


Dependencies:
    bwa
    samtools
    bedtools
"""


import argparse
import csv
import os
import sys
import subprocess


def run_command(cmd):

    print("\nRunning:")
    print(cmd)

    result = subprocess.run(
        cmd,
        shell=True,
        executable="/bin/bash"
    )

    if result.returncode != 0:
        sys.exit(f"\nERROR: command failed:\n{cmd}")


def calculate_copy_number(depth_file, output, query, subject, sample):

    query_depth = None
    query_name = None
    subjects = []


    with open(depth_file) as f:

        reader = csv.reader(
            f,
            delimiter="\t"
        )

        for row in reader:

            if len(row) < 5:
                continue

            contig = row[0]
            annotation = row[3]

            try:
                depth = float(row[4])
            except ValueError:
                continue


            annotation_lower = annotation.lower()


            # Find query gene (first match)
            if (
                query_depth is None
                and query.lower() in annotation_lower
            ):
                query_depth = depth
                query_name = query


            # Find all subject genes
            if subject.lower() in annotation_lower:

                subjects.append(
                    (contig, depth)
                )


    if query_depth is None:
        sys.exit(
            f"ERROR: Query gene '{query}' not found."
        )


    if len(subjects) == 0:
        sys.exit(
            f"ERROR: Subject gene '{subject}' not found."
        )



    with open(output, "w", newline="") as out:

        writer = csv.writer(
            out,
            delimiter="\t"
        )


        writer.writerow([
            "sample",
            "query_gene",
            "query_depth",
            "subject_contig",
            "subject_gene",
            "subject_depth",
            "copy_number"
        ])



        for contig, subject_depth in subjects:

            copy_number = subject_depth / query_depth


            writer.writerow([

                sample,

                query_name,

                f"{query_depth:.2f}",

                contig,

                subject,

                f"{subject_depth:.2f}",

                f"{copy_number:.3f}"

            ])



def main():

    parser = argparse.ArgumentParser(
        description="Calculate gene copy number from read depth"
    )


    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help="Genome FASTA file"
    )


    parser.add_argument(
        "-1",
        "--reads1",
        required=True,
        help="Forward reads"
    )


    parser.add_argument(
        "-2",
        "--reads2",
        required=True,
        help="Reverse reads"
    )


    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        help="Genome annotation GFF"
    )


    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV"
    )


    parser.add_argument(
        "-q",
        "--query",
        required=True,
        help="Reference single-copy gene (e.g. rpoB)"
    )


    parser.add_argument(
        "-s",
        "--subject",
        required=True,
        help="Target gene (e.g. OXA-23)"
    )


    args = parser.parse_args()



    fasta = args.fasta


    # Extract sample name from sample.fa
    sample = os.path.splitext(
        os.path.basename(fasta)
    )[0]


    bam = sample + ".bam"
    bed = sample + ".genes.bed"
    depth = sample + ".gene_depth.tsv"



    # Check input files

    for file in [
        fasta,
        args.reads1,
        args.reads2,
        args.gff
    ]:

        if not os.path.exists(file):

            sys.exit(
                f"ERROR: Cannot find file: {file}"
            )



    # --------------------------------------------------
    # 1. Index genome
    # --------------------------------------------------

    run_command(
        f"bwa index {fasta}"
    )



    # --------------------------------------------------
    # 2. Map reads
    # --------------------------------------------------

    run_command(
        f"bwa mem -t 16 {fasta} "
        f"{args.reads1} {args.reads2} | "
        f"samtools sort -o {bam}"
    )


    run_command(
        f"samtools index {bam}"
    )



    # --------------------------------------------------
    # 3. Convert CDS features from GFF to BED
    # --------------------------------------------------

    run_command(
        f"""awk -F'\\t' 'BEGIN{{OFS="\\t"}}
        $3=="CDS" {{
            gsub(/ /,"_",$9);
            print $1,$4-1,$5,$9
        }}' {args.gff} > {bed}"""
    )



    # --------------------------------------------------
    # 4. Calculate mean depth per gene
    # --------------------------------------------------

    run_command(
        f"bedtools coverage "
        f"-a {bed} "
        f"-b {bam} "
        f"-mean > {depth}"
    )



    # --------------------------------------------------
    # 5. Calculate copy number
    # --------------------------------------------------

    calculate_copy_number(
        depth_file=depth,
        output=args.output,
        query=args.query,
        subject=args.subject,
        sample=sample
    )



    print("\nFinished!")
    print(f"Sample: {sample}")
    print(f"Output: {args.output}")



if __name__ == "__main__":
    main()
