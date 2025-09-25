import sys
import dataclasses
import subprocess
import gzip
from Bio import SeqIO
from typing import Generator


@dataclasses.dataclass(frozen=True)
class ReadnameVarTags:
    readname: str
    vartag: str
    remain_fields: str

def find_vartag_readnames(lookfor_vartag: str, readname_file: str) -> set[str]:
    """
     Find read names with a specific variant tag.
 
     Args:
     - lookfor_vartag (str): The variant tag to look for.
     - readname_file (str): Path to the file containing read names and variant tags.
 
     Returns:
     - set[str]: A set of read names that match the specified variant tag.
     """
    matched_readnames = set()
    with open(readname_file, 'rt') as fh:
        for line in fh:
            rn, vartag, *_ = line.split('\t')
            if vartag and vartag.startswith(lookfor_vartag):
                matched_readnames.add(rn)
    return matched_readnames


def extract_variant_seqs(readname_set: set[str], fastq_file: str) -> dict[str, str]:
    """
    Extracts sequences from a FASTQ file for a given set of read names.

    Args:
        readname_set (set): A set of read names to extract sequences for.
        fastq_file (str): The path to the FASTQ file. The file can be gzipped.

    Returns:
        Dict[str, str]: A dictionary with read names as keys and sequences as values.
    """
    if fastq_file.endswith('.gz'):
        with gzip.open(fastq_file, "rt") as fq:
            return dict((rec.id, rec.seq) for rec in SeqIO.parse(fq, 'fastq') if rec.id in readname_set)

    with open(fastq_file, "rt") as fq:
        return dict((rec.id, rec.seq) for rec in SeqIO.parse(fq, 'fastq') if rec.id in readname_set)


def write_sequence_dict_to_fasta(matched_fq: dict[str, str], out_fa_name: str) -> None:
    """
    Writes sequences from a dictionary to a FASTA file.

    Args:
        matched_fq (Dict[str, str]): Dictionary with read names as keys and sequences as values.
        out_fa_name (str): Output FASTA file name.
    """
    with open(out_fa_name, 'w') as ofa:
        for identifier, seq in matched_fq.items():
            fasta_record = f">{identifier}\n{seq}\n"
            ofa.write(fasta_record)


def create_ref_fasta(reference_seq: str, reference_identifier: str) -> None:
    """
    Creates a reference FASTA file with the given sequence and identifier.

    Args:
        reference_seq (str): The reference sequence to write to the FASTA file.
        reference_identifier (str): The identifier for the reference sequence.
    """
    with open(f'{reference_identifier}.fa', 'w') as fh:
        fh.write(f'>{reference_identifier}\n{reference_seq}\n')


def align_sort(reference_fasta: str, input_fasta: str, sam_file_name: str) -> None:
    """
    Aligns and sorts sequences using minimap2 and samtools.

    Args:
        reference_fasta (str): The path to the reference FASTA file.
        input_fasta (str): The path to the input FASTA file.
        sam_file_name (str): The base name for the output SAM file.
    """
    align_sort_cmd = f'minimap2 -ax map-hifi -t 8 {reference_fasta} {input_fasta} | samtools sort -o {sam_file_name}.sam'
    process = subprocess.run(align_sort_cmd, shell=True)
    if process.returncode == 0:
        print(
            f'Alignment and sorting completed successfully. {sam_file_name}.sam returned')
    else:
        print(f'Failed on alignment and sorting: {align_sort_cmd}.')
        sys.exit(process.returncode)


def parse_readname_record(line: str) -> ReadnameVarTags:
    readname, vartag, *others = line.split('\t')
    return ReadnameVarTags(readname, vartag, others)


def should_update_vartag(record: ReadnameVarTags, original_vartag: str, true_reads: set[str], invert_condition: bool) -> bool:
    return (
        record.vartag
        and record.vartag.startswith(original_vartag)
        and (record.readname in true_reads) != invert_condition)


def modify_vartags(readname_file: str, 
                   true_reads: set[str],
                   original_vartag: str = 'HDR', 
                   replacement_vartag: str = 'Non-HDR-without-ITR',
                   invert_condition: bool = True) -> Generator[str, str, str]:
    """
    Updates variant tags in the readname file based on true HDR reads.

    Args:
        readname_file (str): Path to the readname file.
        true_reads (set[str]): Set of true vartag read names.
        original_vartag (str): Original variant tag to replace.
        replacement_vartag (str): Replacement variant tag.
        invert_condition (bool): If True, the condition will be inverted.

    Yields:
        Generator[str, str, str]: Updated readname records with modified variant tags.
    """
    with open(readname_file, 'r') as fh:
        lines = fh.readlines()
        for line in lines:
            record = parse_readname_record(line)
    
            if should_update_vartag(record, original_vartag, true_reads, invert_condition):
                updated_vartag = record.vartag.replace(original_vartag, replacement_vartag)
                yield dataclasses.replace(record, vartag=updated_vartag)
            else:
                yield record

def write_output(readname_record: Generator[str, str, str], filename:str) -> None:
    """
    Writes updated readname records to an output file.

    Args:
        readname_record (Generator[str, str, str]): Generator yielding updated readname records.
        filename (str): filename for readname txt file.
    """
    with open(f"readname_{filename}.txt", 'w') as fh:
        for rec in readname_record:
            readname = rec.readname
            vartag = rec.vartag
            remained = rec.remain_fields
            fh.write('\t'.join([readname, vartag, *remained]))


