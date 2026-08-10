#!/usr/bin/env python3

"""
Usage:

python gene_copy_number.py \
    -f sample.fa \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -g sample.gff \
    -o output.tsv \
    -q <gene_name_of_reference> \
    -s <gene_name_of_interest>
    -t <number_of_threads>

Also accepts GFF3 input (e.g. -g sample.gff3), including files with an
embedded ##FASTA sequence block at the end (as produced by Bakta/Prokka).

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

    # "set -o pipefail" ensures that if any command in a pipe fails
    # (e.g. bwa mem), the whole pipeline is reported as failed, even
    # if a downstream command (e.g. samtools sort) succeeds on empty
    # input. Without this, a failed mapping step can silently produce
    # a valid-but-empty BAM, which later shows up as all-zero depth.
    result = subprocess.run(
        "set -o pipefail; " + cmd,
        shell=True,
        executable="/bin/bash"
    )

    if result.returncode != 0:
        sys.exit(f"\nERROR: command failed:\n{cmd}")


def get_bam_contigs(bam):

    result = subprocess.run(
        f"samtools view -H {bam}",
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        sys.exit(f"\nERROR: could not read BAM header from {bam}")

    contigs = set()

    for line in result.stdout.splitlines():

        if line.startswith("@SQ"):

            for field in line.split("\t"):

                if field.startswith("SN:"):
                    contigs.add(field[3:])

    return contigs


def get_bed_contigs(bed):

    contigs = set()

    with open(bed) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            contigs.add(line.split("\t")[0])

    return contigs


def validate_contigs_match(bam, bed):

    bam_contigs = get_bam_contigs(bam)
    bed_contigs = get_bed_contigs(bed)

    overlap = bam_contigs & bed_contigs

    if len(overlap) == 0:

        sys.exit(
            "\nERROR: No contig names in common between the BAM "
            f"({bam}) and the annotation BED ({bed}).\n"
            "This means the genome FASTA used for mapping (-f) and "
            "the GFF/GFF3 used for annotation (-g) do not match - "
            "coverage will silently come out as all zeros otherwise.\n\n"
            f"Example BAM contig(s): {sorted(bam_contigs)[:3]}\n"
            f"Example BED contig(s): {sorted(bed_contigs)[:3]}\n\n"
            "Fix: make sure -f is the exact same assembly (same contig "
            "names) as the one that was annotated to produce -g. If your "
            "GFF3 has an embedded ##FASTA section (common with Bakta/"
            "Prokka), extract and use that as -f instead, e.g.:\n"
            "  awk '/^##FASTA/{found=1; next} found' your.gff3 > "
            "genome_from_gff3.fa"
        )

    elif len(overlap) < len(bed_contigs):

        print(
            "\nWARNING: only "
            f"{len(overlap)}/{len(bed_contigs)} annotation contigs "
            "were found in the BAM. Depth for genes on the missing "
            "contigs will be reported as 0."
        )


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
        help="Genome annotation GFF or GFF3 file"
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


    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="Number of threads for bwa mem (default: 16)"
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
        f"bwa mem -t {args.threads} {fasta} "
        f"{args.reads1} {args.reads2} | "
        f"samtools sort -o {bam}"
    )


    run_command(
        f"samtools index {bam}"
    )



    # --------------------------------------------------
    # 3. Convert CDS features from GFF/GFF3 to BED
    #
    # Works for both plain GFF and GFF3:
    #   - GFF3 attribute strings (key=value;key=value) are kept as-is
    #     in column 9, which is what the query/subject substring
    #     matching in calculate_copy_number() relies on.
    #   - Lines starting with "#" (GFF3 pragmas/comments, e.g.
    #     ##gff-version 3, ##sequence-region) are skipped explicitly.
    #   - Processing stops at a "##FASTA" line, so any embedded
    #     genome sequence appended after the annotations (common in
    #     Bakta/Prokka GFF3 output) is never scanned.
    # --------------------------------------------------

    run_command(
        f"""awk -F'\\t' 'BEGIN{{OFS="\\t"}}
        /^##FASTA/ {{exit}}
        /^#/ {{next}}
        $3=="CDS" {{
            gsub(/ /,"_",$9);
            print $1,$4-1,$5,$9
        }}' {args.gff} > {bed}"""
    )



    # --------------------------------------------------
    # 3.5. Sanity check: do the BAM and BED contig names
    #      actually overlap? If not, bedtools coverage
    #      would silently report 0 for every gene.
    # --------------------------------------------------

    validate_contigs_match(bam, bed)



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
