#!/usr/bin/env python3
# coding: utf-8

# import modules
import os
import sys
import gzip
from typing import Set, List

from Bio import SeqIO


def read_names_to_set(txt_file: str) -> List:
    """Read the list of read names from a text file"""
    with open(txt_file, 'rt') as readname_file:
        return set(line.strip() for line in readname_file)


def extract_reads_by_readnames(read_names: Set, input_fastq: str) -> None:
    """Iterate through the FASTQ.gz file matching reads to the output"""
    filename = os.path.basename(input_fastq).split('.')[0]

    with gzip.open(input_fastq, 'rt') as i_fq:
        matching_record = []
        for rec in SeqIO.parse(i_fq, "fastq"):
            if rec.id in read_names:
                rec.id = f"{rec.id}_{filename}"
                rec.description = ''
                rec.name = f"{rec.id}_{filename}"
                matching_record.append(rec)
        return matching_record


def write_records_to_fq(seq_record: List, output_fq: str) -> None:
    with open(output_fq, "wt") as ofq:
        for rec in seq_record:
            SeqIO.write(rec, ofq, "fastq")
    return 0


if __name__ == "__main__":
    read_names_file = sys.argv[1]
    input_fastq_file = sys.argv[2]
    output_fastq_file = sys.argv[3]
    sample_name = os.path.basename(input_fastq_file)
    print(f'>{"Start to process: " + sample_name :=^120}<')

    read_names = read_names_to_set(read_names_file)
    matching_records = extract_reads_by_readnames(read_names, input_fastq_file)
    write_records_to_fq(matching_records, output_fastq_file)

    print(
        f"Extracted and Modified {len(matching_records)} reads and saved them in {output_fastq_file}.")
