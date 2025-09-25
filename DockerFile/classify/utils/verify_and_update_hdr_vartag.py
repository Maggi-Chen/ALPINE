import pysam
from utils import data_utils as du

def find_true_hdr_reads(sam_file: str, ha_length: int, threshold: float = 0.96) -> set[str]:
    """
    Identifies HDR reads from a SAM file based on alignment length.

    Args:
        sam_file (str): Path to the SAM file.
        ha_length (int): Length of the homology arm (HA) sequence.
        threshold (float): Minimum alignment ratio to consider a read as HDR (default is 0.96).

    Returns:
        set[str]: Set of read names that qualify as HDR.
    """
    hdr_readnames = set()
    edit_distance = round(ha_length * threshold, 0)

    with pysam.AlignmentFile(sam_file, 'r') as sam:
        for read in sam:
            if read.is_secondary or read.has_tag("SA") or read.is_supplementary:
                continue

            aligned_pairs = read.get_aligned_pairs(matches_only=True)
            aligned_length = len(aligned_pairs)

            if aligned_length >= edit_distance:
                hdr_readnames.add(read.query_name)

    return hdr_readnames


def verify_update_hdr_tag(
    readname_tsv: str, 
    filtered_fastq: str,   
    left_ha_seq: str, 
    right_ha_seq: str, 
    output_filename: str,
    lookfor_vartag: str = 'HDR', 
    ha_match_threshold: float = 0.96
) -> None:
    """
    Extract variant sequences from the filtered FASTQ file based on the read names
    that contain the specified variant tag (e.g., HDR). These sequences are then
    written to a FASTA file for further alignment and analysis.
    """
    # extract variant sequences from fastq and output those sequences into fasta
    unfiltered_hdr_readname_set = du.find_vartag_readnames(
        lookfor_vartag, readname_tsv)
    variant_seq_dict = du.extract_variant_seqs(
        unfiltered_hdr_readname_set, filtered_fastq)
    du.write_sequence_dict_to_fasta(variant_seq_dict, f'{lookfor_vartag}.fa')

    # align & sort
    ha_names = ('left_ha', 'right_ha')
    ha_seqs = (left_ha_seq, right_ha_seq)

    readname_dict = {}
    # Process each HA sequence
    for ha_name, ha_seq in zip(ha_names, ha_seqs):
        # Create reference FASTA for the HA sequence
        du.create_ref_fasta(ha_seq, ha_name)

        # Align the HA ref with the variant read fasta
        du.align_sort(f'{ha_name}.fa', f'{lookfor_vartag}.fa', f'{lookfor_vartag}_{ha_name}')

        # Extract readnames
        ha_readnames = find_true_hdr_reads(f'{lookfor_vartag}_{ha_name}.sam', len(ha_seq), threshold=ha_match_threshold)

        # Store the readnames for this HA
        readname_dict[ha_name] = ha_readnames

    # Find common readnames in both left and right HAs
    true_hdr_readnames = readname_dict['left_ha'] & readname_dict['right_ha']

    # Update HDR variant tag for the readname file
    modified_lines = du.modify_vartags(readname_tsv, true_hdr_readnames)
    du.write_output(modified_lines, output_filename)

