import subprocess
from utils import data_utils as du
import sys

def make_blast_database(db_name: str) -> None:
    """
    Create a BLAST database from a sequence.

    Args:
        itr_seq (str): The sequence to create the database from.
        db_name (str): The name of the database.
    """
    make_db_cmd = f'makeblastdb -in {db_name}.fa -dbtype nucl -parse_seqids -out {db_name}.db'
    
    process = subprocess.run(make_db_cmd, shell=True)
    if process.returncode == 0:
        print(f'BLAST database created successfully: {db_name}.db')
    else:
        print(f'Failed to create BLAST database: {make_db_cmd}')
        sys.exit(process.returncode)


def extract_itr_readnames(fasta_file: str, db_name: str, seed_size: int = 15, perc_identity: int = 70) -> set[str]:
    """
    Perform blast, then extract read names from BLAST stdout.

    Args:
        fasta_file (str): Path to the FASTA file.
        db_name (str): The name of the BLAST database.
        seed_size (int): The word size for the BLAST search (default is 16).
        perc_identity (int): The minimum percentage identity for a hit (default is 50%).

    Returns:
        set[str]: Set of read names extracted from the BLAST output.
    """
    readnames = set()

    # 6 is tabular output for hits
    blast_cmd = f'blastn -num_threads 15 -word_size {seed_size} -perc_identity {perc_identity} -query {fasta_file} -db {db_name}.db -outfmt 6'
    blast_output = subprocess.run(
        blast_cmd, shell=True, capture_output=True, text=True)
    if blast_output.returncode == 0:
        print(f'BLAST completed successfully.')
        for line in blast_output.stdout.splitlines():
            readnames.add(line.split('\t')[0])
        return readnames      
    else:
        print(f'Failed on blast {blast_cmd}.')
        sys.exit(blast_output.returncode)


def verify_itr_seq(readname_tsv: str,
                    filtered_fastq: str,
                    output_filename: str,
                    left_itr_seq: str,
                    right_itr_seq: str,
                    lookfor_vartag: str = 'Non-HDR-without-ITR',
                    seed_size: int = 15,
                    percent_identity: int = 70) -> None:
    """
    Verify and update the variant tag for reads based on the presence of ITR sequences in read.

    Args:
        readname_tsv (str): Path to the TSV file containing read names.
        filtered_fastq (str): Path to the filtered FASTQ file.
        output_filename (str): Path to the output file where results will be written.
        left_itr_seq (str): Sequence of the left ITR.
        right_itr_seq (str): Sequence of the right ITR.
        lookfor_vartag (str): Variant tag to look for in the read names (default is 'Non-HDR-without-ITR').

    """
    
    # extract variant sequences from fastq and output those sequences into fasta
    unfiltered_itr_readname_set = du.find_vartag_readnames(
        lookfor_vartag, readname_tsv)
    
    variant_seq_dict = du.extract_variant_seqs(
        unfiltered_itr_readname_set, filtered_fastq)
    
    du.write_sequence_dict_to_fasta(variant_seq_dict, f'{lookfor_vartag}.fa')

    # Loop each ITR sequence to create a BLAST database, blast, and extract found readnames
    itr_names = ('left_itr', 'right_itr')
    itr_seqs = (left_itr_seq, right_itr_seq)

    readnames_with_itr = set()

    # Process each side ITR sequence
    for itr_name, itr_seq in zip(itr_names, itr_seqs):

        du.create_ref_fasta(itr_seq, itr_name)

        make_blast_database(itr_name)

        current_readname_set = extract_itr_readnames(
            f'{lookfor_vartag}.fa', itr_name, seed_size=seed_size, perc_identity=percent_identity)
        readnames_with_itr.update(current_readname_set)

   # Process each line in readname file and correct the vartag if the readname is in readnames_with_itr 
    modified_lines = du.modify_vartags(readname_tsv,
                                       readnames_with_itr,
                                       'Non-HDR-without-ITR',
                                       'Non-HDR-with-ITR',
                                       invert_condition=False)
    
    du.write_output(modified_lines, output_filename)
