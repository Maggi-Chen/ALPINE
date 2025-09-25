#!/usr/bin/env python3
# coding: utf-8

# import modules
import argparse
import gzip
import os
import time
from datetime import datetime
from typing import NamedTuple
from utils.create_size_quality_plots import create_combined_histograms
from Bio import SeqIO

class Args(NamedTuple):
    """ Command-line arguments"""
    input: str
    output:str
    forw: str
    reve: str
    min_qual: int

def get_args() -> Args:
    """ Get command-line arguments"""

    parser = argparse.ArgumentParser(usage='filtering long-reads fq.gz file.',
                                    description="filtering long-reads fq.gz file based \
                                        read length and average read base quality score.\
                                        Output file are filtered and saved as gz format.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                    )
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="input fq.gz file")
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="output directory", default="")
    parser.add_argument("--sample", type=str, required=False,
                        help="sample name used for naming output file", default="Sample")
    parser.add_argument("-f", "--forw", type=str, required=False, default=None,
                        help="forward primer sequence")
    parser.add_argument("-r", "--reve", type=str, required=False, default=None,
                        help="reverse primer sequence")
    parser.add_argument("-q", "--min_qual", type=float, default=30,
                        help="minimal average read base quality")
    args = parser.parse_args()
    return args

def main() -> None:
    args = get_args()

    # Get current date and time
    current_datetime = datetime.now()
    current_date = current_datetime.date()
    current_time_formatted = current_datetime.strftime("%H:%M:%S")

    print(
        f"Program start: {current_date}, {current_time_formatted}.")
    
    start_time = time.time()
    (readlen_orig, readlen_filtered, readqual_orig, readqual_filtered) = filter_fastqgz(args.input,
                   args.output,
                   args.sample,
                args.forw,
                revcom(args.reve),
                args.min_qual)


    create_combined_histograms(readlen_orig, readqual_orig, args.sample+"_original")
    create_combined_histograms(readlen_filtered, readqual_filtered, args.sample+"_filtered")
    print(f'Run complete. Total elapsed time:{(time.time() - start_time)/60:.1f} \
            minutes.')
    
def avgqual(readqual:SeqIO) -> int:
    """Calculate average phred score in a read.
    Args:
        readqual (SeqIO): read quality object from SeqIO sequece record
    Returns:
        either 0 or mean of quality score
    """
    numbases = len(readqual)
    if numbases < 100:
        return 0
    else:
        return int(sum(readqual) / numbases)

def revcom(seq:str) -> str:
    """Reverse, complement a dna sequence
    """
    if not seq:
        return seq
    comp = dict(zip(list('ATGC'), list('TACG')))
    return ''.join(
        base if base not in comp else comp[base] for base in reversed(seq))

def is_subsequence(seq, subseq):
    '''Return True if 10bp of any subseq matches with seq'''
    if not subseq:
        return True
    for i in range(len(subseq)-9):
        targetseq=subseq[i:i+10]
        if targetseq in seq:
            return True
    return False


def primer_in(readseq:str, primerf:str, primerr:str, revprimerf:str, revprimerr:str) -> bool:
    """Test if first or last 20 base nucleotides matched to either forward or reverse primers.
    """
    head = readseq[:50]
    tail = readseq[-50:]

    if is_subsequence(head,primerf) and is_subsequence(tail,primerr):
        return True
    if is_subsequence(head,revprimerr) and  is_subsequence(tail,revprimerf):
        return True
    return False


def filter_fastqgz(fastqpath:str, output_dir:str, sample:str, primer_f:str, primer_r:str, qual=30) -> None:
    """
    Filter a FASTQ file based on quality, forward, and reverse primer matches.
    Args:
        fastqpath (str): Path to the input FASTQ file in gzip-compressed format.
        output_dir (str): Directory where filtered FASTQ file and log will be saved.
        primer_f (str): Forward primer sequence.
        primer_r (str): Reverse primer sequence.
        qual (int): Minimum quality score threshold. 
                    Reads with lower average quality will be discarded.

    Returns:
        0 if the filtering process completes.
    """

    input_file_name = os.path.basename(fastqpath).split('.fastq')[0]
    output_file_name = sample+'.fastq.gz'
    output_fq_path = os.path.join(output_dir, output_file_name)
    rev_primer_f = revcom(primer_f)
    rev_primer_r = revcom(primer_r)
    cut_off_qual = qual if qual else 30

    low_qual = 0
    fulllen = 0
    nonfulllen = 0

    readlen_orig = []; readlen_filtered = []
    readqual_orig = []; readqual_filtered = []
   
    try:
        if fastqpath.endswith('.gz'):
            in_f = gzip.open(fastqpath, 'rt')
        else:
            in_f = open(fastqpath, 'r')
        with gzip.open(output_fq_path, 'wt') as out_f:
            for rec in SeqIO.parse(in_f, "fastq"):
                seq = rec.seq
                phred_score = rec.letter_annotations["phred_quality"]

                readqual_orig += [avgqual(phred_score)]
                readlen_orig += [len(seq)]
                if avgqual(phred_score) < cut_off_qual:
                    low_qual += 1
                    continue
                readqual_filtered += [avgqual(phred_score)]

                if primer_in(seq,
                             primer_f,
                             primer_r,
                             rev_primer_f,
                             rev_primer_r):
                    SeqIO.write(rec, out_f, "fastq")
                    fulllen += 1
                    readlen_filtered += [len(seq)]
                    if fulllen % 1000 == 0: 
                        print(f"number of reads being filtered: {fulllen}")

                else:
                    nonfulllen += 1
    

        log_file_name = os.path.join(output_dir,
                                        sample+".filtering.log")
        with open(log_file_name, 'a') as lf:
            lf.write('Filter Read: ' + input_file_name
                        + '\tFulllen Read (pass filter): ' + str(fulllen)
                        + '\tNon-fulllen Read (No primer): ' + str(nonfulllen)
                        + '\tLow-quality Read: ' + str(low_qual) + '\n')

    except Exception as e:
        print(f"An error occurred: {str(e)}")
    
    return (readlen_orig, readlen_filtered, readqual_orig, readqual_filtered)

if __name__ == "__main__":
    main()

