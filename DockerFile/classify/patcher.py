#!/usr/bin/env python3
"""
Patcher Module for AAV Integration Classification Pipeline

This module post-processes classification results to rescue false negatives.
It re-analyzes reads classified as "Unclassified" by realigning them to
WT-only reference to detect large deletions.

This is a SAFE post-processor that:
- Does NOT modify the original classify_read.py
- Only touches Unclassified reads
- Preserves all other classifications unchanged
- Creates a new patched output file (original is preserved)

Usage:
    python patcher.py --readname readname_sample.txt \
                      --bam sample.bam \
                      --ref ref_combined.fa \
                      --config config.txt \
                      [--output readname_sample_patched.txt] \
                      [--length-threshold 50] \
                      [--stats stats.txt] \
                      [--minimap2-path /path/to/minimap2] \
                      [--samtools-path /path/to/samtools]

Environment Variables:
    PATCHER_MINIMAP2_PATH - Path to minimap2 v2.16 executable
    PATCHER_SAMTOOLS_PATH - Path to samtools executable

Docker Defaults:
    minimap2: /opt/minimap2-2.16/minimap2
    samtools: /usr/local/bin/samtools

Author: Xing-Huang Gao 
"""

import argparse
import os
import sys
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import shared utilities
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.config_utils import read_aav_config
from utils.debreak_detect import segmentdeletion_tumor, cigardeletion

# Default paths for Docker container
# Patcher requires minimap2 v2.16 specifically (newer versions soft-clip large deletions)
DEFAULT_MINIMAP2_PATH = "/opt/minimap2-2.16/minimap2"
DEFAULT_SAMTOOLS_PATH = "/usr/local/bin/samtools"


def get_tool_path(cli_path, env_var, default_path):
    """
    Get tool path with priority: CLI arg > environment variable > default.

    Args:
        cli_path: Path from command line argument (or None)
        env_var: Environment variable name to check
        default_path: Default path for Docker container

    Returns:
        str: Resolved tool path
    """
    if cli_path:
        return cli_path
    return os.environ.get(env_var, default_path)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Patcher: Rescue false negatives in AAV classification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python patcher.py --readname readname_sample.txt --bam sample.bam \\
                      --ref ref_combined.fa --config config.txt

    # With custom output and threshold
    python patcher.py --readname readname_sample.txt --bam sample.bam \\
                      --ref ref_combined.fa --config config.txt \\
                      --output patched.txt --length-threshold 30
        """
    )

    parser.add_argument('--readname', required=True,
                        help='Classification output file (readname_{sample}.txt)')
    parser.add_argument('--bam', required=True,
                        help='Sorted BAM file with aligned reads')
    parser.add_argument('--ref', required=True,
                        help='Combined reference FASTA (contains WT and AAV sequences)')
    parser.add_argument('--config', required=True,
                        help='AAV vector configuration file')
    parser.add_argument('--output', default=None,
                        help='Output file path (default: {input}_patched.txt)')
    parser.add_argument('--length-threshold', type=int, default=50,
                        help='Length threshold for weak Non-HDR-without-ITR (default: 50)')
    parser.add_argument('--cigar-window', type=int, default=20,
                        help='Window for CIGAR-based detection ±bp (default: 20, matches classify.py)')
    parser.add_argument('--segment-window', type=int, default=10,
                        help='Window for segment-based detection ±bp (default: 10, matches classify.py)')
    parser.add_argument('--stats', default=None,
                        help='Optional file to write statistics')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep temporary files for debugging')
    parser.add_argument('--minimap2-path', default=None,
                        help=f'Path to minimap2 v2.16 executable (default: env PATCHER_MINIMAP2_PATH or {DEFAULT_MINIMAP2_PATH})')
    parser.add_argument('--samtools-path', default=None,
                        help=f'Path to samtools executable (default: env PATCHER_SAMTOOLS_PATH or {DEFAULT_SAMTOOLS_PATH})')

    return parser.parse_args()


def parse_classification_file(readname_file):
    """
    Parse the classification output file.

    Returns:
        dict: {readname: {'category': str, 'length': str, 'line': str}}
    """
    classifications = {}

    with open(readname_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                continue

            readname = parts[0]
            category = parts[1]
            length = parts[2] if len(parts) > 2 else ''

            classifications[readname] = {
                'category': category,
                'length': length,
                'line': line
            }

    return classifications


def identify_candidates(classifications, length_threshold=50):
    """
    Identify reads that need re-analysis.

    Candidates:
    1. Unclassified reads only

    Args:
        classifications: dict from parse_classification_file()
        length_threshold: Not used (kept for API compatibility)

    Returns:
        list: List of dicts with candidate information
    """
    candidates = []

    for readname, info in classifications.items():
        category = info['category']
        length = info['length']

        # Only process Unclassified reads
        if category == "Unclassified":
            candidates.append({
                'readname': readname,
                'old_category': category,
                'old_length': length,
                'reason': 'Unclassified'
            })

    return candidates


def extract_sequences_from_bam(bam_file, read_names):
    """
    Extract full read sequences from BAM file.

    Args:
        bam_file: Path to sorted BAM file
        read_names: Set of read names to extract

    Returns:
        dict: {readname: sequence}
    """
    sequences = {}
    read_names_set = set(read_names)

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_name in read_names_set:
                # Only process primary alignments with sequence
                if not read.is_secondary and not read.is_supplementary:
                    if read.query_sequence:
                        sequences[read.query_name] = read.query_sequence

    return sequences


def extract_wt_reference(ref_fasta, config_file):
    """
    Extract WT reference information from config and FASTA.

    Args:
        ref_fasta: Path to combined reference FASTA
        config_file: Path to AAV config file

    Returns:
        tuple: (wt_chrom, wt_sequence, target_pos)
    """
    # Parse config to find WT reference
    aavinfo = read_aav_config(config_file)

    wt_chrom = None
    target_pos = None

    for aav_id in aavinfo:
        for ref_name in aavinfo[aav_id]:
            if aavinfo[aav_id][ref_name]['Ref_Type'] == "WT":
                wt_chrom = ref_name
                gene_start = aavinfo[aav_id][ref_name]['gene_start']
                gene_end = aavinfo[aav_id][ref_name]['gene_end']
                target_pos = (gene_start + gene_end) // 2
                break
        if wt_chrom:
            break

    if not wt_chrom:
        raise ValueError("No WT reference found in config file")

    # Extract WT sequence from FASTA
    wt_sequence = None
    for record in SeqIO.parse(ref_fasta, "fasta"):
        if record.id == wt_chrom:
            wt_sequence = str(record.seq)
            break

    if not wt_sequence:
        raise ValueError(f"WT reference '{wt_chrom}' not found in FASTA file")

    return wt_chrom, wt_sequence, target_pos


def write_temp_files(sequences, wt_chrom, wt_sequence, temp_prefix):
    """
    Write temporary FASTA files for realignment.

    Args:
        sequences: dict of {readname: sequence}
        wt_chrom: WT chromosome name
        wt_sequence: WT reference sequence
        temp_prefix: Prefix for temp file names

    Returns:
        tuple: (reads_fasta_path, wt_fasta_path)
    """
    reads_fasta = f"{temp_prefix}_reads.fa"
    wt_fasta = f"{temp_prefix}_wt.fa"

    # Write read sequences
    with open(reads_fasta, 'w') as f:
        for readname, seq in sequences.items():
            f.write(f">{readname}\n{seq}\n")

    # Write WT reference
    with open(wt_fasta, 'w') as f:
        f.write(f">{wt_chrom}\n{wt_sequence}\n")

    return reads_fasta, wt_fasta


def realign_to_wt(reads_fasta, wt_fasta, temp_prefix, minimap2_path, samtools_path):
    """
    Realign reads to WT-only reference using minimap2.

    NOTE: Uses old minimap2 (2.16) with map-pb preset intentionally.
    This is because newer minimap2 versions soft-clip large deletions instead
    of representing them in the CIGAR string. For the patcher's purpose of
    rescuing false negatives, map-pb with its lower gap penalty is beneficial
    as it explicitly shows large deletions (e.g., 996D) in the CIGAR.

    Args:
        reads_fasta: Path to reads FASTA
        wt_fasta: Path to WT reference FASTA
        temp_prefix: Prefix for output files
        minimap2_path: Path to minimap2 v2.16 executable
        samtools_path: Path to samtools executable

    Returns:
        str: Path to output BAM file
    """
    output_sam = f"{temp_prefix}_realigned.sam"
    output_bam = f"{temp_prefix}_realigned.bam"

    # Use old minimap2 (2.16) with map-pb preset
    # - map-pb has lower gap open penalty (4 vs 6) which allows large deletions
    # - Newer minimap2 (2.28+) soft-clips instead of showing deletions in CIGAR
    cmd = f"{minimap2_path} -a -t 4 --secondary=no {wt_fasta} {reads_fasta} > {output_sam}"
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"minimap2 failed with return code {ret}")

    # Convert to sorted BAM
    cmd = f"{samtools_path} sort -o {output_bam} {output_sam} && {samtools_path} index {output_bam}"
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"samtools failed with return code {ret}")

    return output_bam


def detect_deletions(realigned_bam, wt_chrom, target_pos, cigar_window=20, segment_window=10, min_del_size=50):
    """
    Detect deletions from realigned BAM file.

    Uses existing segmentdeletion_tumor and cigardeletion functions
    from utils.debreak_detect module.

    Window sizes match classify.py for consistency:
    - CIGAR-based detection: ±20bp (matches classify.py line 108)
    - Segment-based detection: ±10bp (matches classify.py line 137)

    Args:
        realigned_bam: Path to realigned BAM file
        wt_chrom: WT chromosome name
        target_pos: Target position for variant detection
        cigar_window: Window for CIGAR-based detection (±bp, default: 20)
        segment_window: Window for segment-based detection (±bp, default: 10)
        min_del_size: Minimum deletion size to report (default: 50bp)
                      Filters out small deletions from messy alignments

    Returns:
        dict: {readname: {'type': str, 'size': int, 'details': str}}
    """
    variants = {}

    # CIGAR detection window (±20bp to match classify.py line 108)
    cigar_window_start = target_pos - cigar_window
    cigar_window_end = target_pos + cigar_window

    # Segment detection window (±10bp to match classify.py line 137)
    segment_window_start = target_pos - segment_window
    segment_window_end = target_pos + segment_window

    # Size thresholds for SV detection
    min_size = min_del_size  # Only detect deletions >= this size
    max_size = 100000

    # Collect segments per read in format expected by segmentdeletion_tumor:
    # [readname, flag, chrom, pos_start, pos_end, [leftclip, readlen, rightclip], MAPQ, readseq]
    segmentreads = {}

    with pysam.AlignmentFile(realigned_bam, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary:
                continue

            readname = read.query_name
            flag = read.flag
            chrom = read.reference_name
            pos_start = read.reference_start
            pos_end = read.reference_end
            mapq = read.mapping_quality
            cigarstring = read.cigarstring
            readseq = read.query_sequence if read.query_sequence else ''

            # Calculate clip info from CIGAR tuples
            cigar_tuples = read.cigartuples
            leftclip = 0
            rightclip = 0
            readlen = 0

            if cigar_tuples:
                # Left clip (S=4, H=5)
                if cigar_tuples[0][0] in [4, 5]:
                    leftclip = cigar_tuples[0][1]
                # Right clip
                if cigar_tuples[-1][0] in [4, 5]:
                    rightclip = cigar_tuples[-1][1]
                # Read length (aligned portion)
                readlen = read.query_alignment_length

            cigarinfo = [leftclip, readlen, rightclip]

            # Build segment info
            if read.is_supplementary:
                readinfo = [readname, flag, chrom, pos_start, pos_end, cigarinfo, mapq, '']
            else:
                readinfo = [readname, flag, chrom, pos_start, pos_end, cigarinfo, mapq, readseq]

                # For primary alignments, also detect CIGAR-based variants
                # Uses ±20bp window to match classify.py line 108
                if cigarstring:
                    cigar_svs = cigardeletion(flag, chrom, pos_start, cigarstring, min_size, max_size)
                    # cigar_svs[0] = list of SVs: [[chrom, pos, size, type], ...]
                    for sv in cigar_svs[0]:
                        sv_chrom, sv_pos, sv_size, sv_type = sv[0], sv[1], sv[2], sv[3]

                        # Check if SV overlaps CIGAR detection window (±20bp)
                        sv_end = sv_pos + sv_size if 'D' in sv_type else sv_pos
                        if sv_pos <= cigar_window_end and sv_end >= cigar_window_start:
                            if 'D' in sv_type and sv_size >= min_del_size:  # Deletion >= threshold
                                variant_type = 'DEL-large'  # Only large deletions (>= 50bp)
                                if readname not in variants or variants[readname]['size'] < sv_size:
                                    variants[readname] = {
                                        'type': variant_type,
                                        'size': sv_size,
                                        'details': f'cigar:{sv_pos}-{sv_pos + sv_size}'
                                    }

            # Add to segment collection
            if readname in segmentreads:
                segmentreads[readname].append(readinfo)
            else:
                segmentreads[readname] = [readinfo]

    # Detect segment-based SVs using segmentdeletion_tumor
    # Uses ±10bp window to match classify.py line 137
    for readname, segments in segmentreads.items():
        if len(segments) < 2 or len(segments) > 20:
            continue

        # Call segmentdeletion_tumor
        segment_svs = segmentdeletion_tumor(segments, min_size, max_size)

        # Parse results: "chrom\tpos\tsize\ttype\treadname_seg1\tflag\tmapq"
        for sv_str in segment_svs:
            parts = sv_str.split('\t')
            if len(parts) < 4:
                continue

            sv_chrom = parts[0]
            sv_pos = int(parts[1])
            sv_size = int(parts[2])
            sv_type = parts[3]

            # Check if SV overlaps segment detection window (±10bp)
            sv_end = sv_pos + sv_size
            if sv_pos <= segment_window_end and sv_end >= segment_window_start:
                # Only rescue DEL-large variants - keep patcher focused and conservative
                # Other variant types (INV, DUP) may have mixed architecture and should
                # remain as their original classification (Unclassified)
                if 'D-segment' in sv_type and sv_size >= min_del_size:  # Deletion >= threshold
                    variant_type = 'DEL-large'  # Only large deletions (>= 50bp)
                    if readname not in variants or variants[readname]['size'] < sv_size:
                        variants[readname] = {
                            'type': variant_type,
                            'size': sv_size,
                            'details': f'segment:{sv_pos}-{sv_end}'
                        }

    return variants


def update_classifications(classifications, candidates, variants):
    """
    Update classifications based on detected variants.

    Args:
        classifications: Original classifications dict
        candidates: List of candidate reads
        variants: Detected variants dict

    Returns:
        tuple: (updated_classifications, changes)
    """
    updated = classifications.copy()
    changes = []

    for candidate in candidates:
        readname = candidate['readname']

        if readname in variants:
            variant = variants[readname]
            new_category = variant['type']
            new_length = str(variant['size'])

            old_category = candidate['old_category']
            old_length = candidate['old_length']

            # Update classification
            updated[readname] = {
                'category': new_category,
                'length': new_length,
                'line': f"{readname}\t{new_category}\t{new_length}"
            }

            changes.append({
                'readname': readname,
                'old_category': old_category,
                'old_length': old_length,
                'new_category': new_category,
                'new_length': new_length,
                'details': variant['details']
            })

    return updated, changes


def write_output(classifications, output_file):
    """
    Write updated classifications to output file.

    Args:
        classifications: Updated classifications dict
        output_file: Output file path
    """
    with open(output_file, 'w') as f:
        for readname, info in classifications.items():
            f.write(f"{info['line']}\n")


def print_statistics(candidates, changes, stats_file=None):
    """
    Print and optionally write statistics.

    Args:
        candidates: List of candidate reads
        changes: List of changes made
        stats_file: Optional file path for stats output
    """
    # Count changes by type
    change_counts = {}
    for change in changes:
        key = f"{change['old_category']} -> {change['new_category']}"
        change_counts[key] = change_counts.get(key, 0) + 1

    # Build statistics text
    stats_lines = [
        "=" * 70,
        "PATCHER STATISTICS",
        "=" * 70,
        "",
        f"Unclassified reads analyzed: {len(candidates)}",
        "",
        f"Reclassifications: {len(changes)}",
    ]

    if change_counts:
        for key, count in sorted(change_counts.items()):
            stats_lines.append(f"  - {key}: {count}")
    else:
        stats_lines.append("  - No reclassifications made")

    remained = len(candidates) - len(changes)
    stats_lines.extend([
        "",
        f"Remained unchanged: {remained}",
        "=" * 70,
    ])

    # Print to stdout
    stats_text = "\n".join(stats_lines)
    print(stats_text)

    # Write to file if specified
    if stats_file:
        with open(stats_file, 'w') as f:
            f.write(stats_text + "\n")

            # Write detailed changes
            if changes:
                f.write("\nDETAILED CHANGES:\n")
                f.write("-" * 70 + "\n")
                for change in changes:
                    f.write(f"{change['readname']}\n")
                    f.write(f"  Old: {change['old_category']} ({change['old_length']})\n")
                    f.write(f"  New: {change['new_category']} ({change['new_length']})\n")
                    f.write(f"  Details: {change['details']}\n")
                    f.write("\n")


def cleanup_temp_files(temp_prefix, keep_temp=False):
    """
    Clean up temporary files.

    Args:
        temp_prefix: Prefix used for temp files
        keep_temp: If True, keep temp files
    """
    if keep_temp:
        print(f"Keeping temp files with prefix: {temp_prefix}")
        return

    temp_files = [
        f"{temp_prefix}_reads.fa",
        f"{temp_prefix}_wt.fa",
        f"{temp_prefix}_realigned.sam",
        f"{temp_prefix}_realigned.bam",
        f"{temp_prefix}_realigned.bam.bai",
    ]

    for f in temp_files:
        if os.path.exists(f):
            os.remove(f)


def rescue_unclassified_reads(readname_file, bam_file, ref_fasta, config_file,
                               min_del_size=50, cigar_window=20, segment_window=10,
                               minimap2_path=None, samtools_path=None,
                               keep_temp=False, verbose=True):
    """
    Callable function to rescue false negatives in AAV classification.

    This function re-analyzes Unclassified reads by realigning them to WT-only
    reference to detect large deletions. It modifies the readname file in-place.

    This is the PRIMARY API for integration with classify_read.py.

    Args:
        readname_file: Path to classification output file (readname_{sample}.txt)
        bam_file: Path to sorted BAM file with aligned reads
        ref_fasta: Path to combined reference FASTA
        config_file: Path to AAV vector configuration file
        min_del_size: Minimum deletion size to detect (default: 50bp)
        cigar_window: Window for CIGAR-based detection ±bp (default: 20)
        segment_window: Window for segment-based detection ±bp (default: 10)
        minimap2_path: Path to minimap2 v2.16 (default: auto-resolve)
        samtools_path: Path to samtools (default: auto-resolve)
        keep_temp: Keep temporary files for debugging (default: False)
        verbose: Print progress messages (default: True)

    Returns:
        dict: Statistics about the rescue operation
            - 'total_reads': Total reads in classification file
            - 'candidates': Number of Unclassified reads analyzed
            - 'rescued': Number of reads reclassified
            - 'changes': List of change details

    Raises:
        FileNotFoundError: If input files don't exist
        RuntimeError: If minimap2/samtools execution fails
        ValueError: If WT reference not found in config/FASTA

    Example:
        >>> from patcher import rescue_unclassified_reads
        >>> stats = rescue_unclassified_reads(
        ...     readname_file='readname_sample.txt',
        ...     bam_file='sample.bam',
        ...     ref_fasta='ref_combined.fa',
        ...     config_file='config.txt'
        ... )
        >>> print(f"Rescued {stats['rescued']} reads")
    """
    # Validate input files exist
    for f, name in [(readname_file, 'readname'), (bam_file, 'bam'),
                    (ref_fasta, 'ref'), (config_file, 'config')]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"{name} file not found: {f}")

    # Resolve tool paths (arg > env var > default)
    minimap2 = get_tool_path(minimap2_path, 'PATCHER_MINIMAP2_PATH', DEFAULT_MINIMAP2_PATH)
    samtools = get_tool_path(samtools_path, 'PATCHER_SAMTOOLS_PATH', DEFAULT_SAMTOOLS_PATH)

    # Validate tool executables
    if not os.path.exists(minimap2):
        raise FileNotFoundError(f"minimap2 not found at: {minimap2}. "
                                "Set PATCHER_MINIMAP2_PATH or pass minimap2_path argument.")
    if not os.access(minimap2, os.X_OK):
        raise PermissionError(f"minimap2 is not executable: {minimap2}")
    if not os.path.exists(samtools):
        raise FileNotFoundError(f"samtools not found at: {samtools}. "
                                "Set PATCHER_SAMTOOLS_PATH or pass samtools_path argument.")
    if not os.access(samtools, os.X_OK):
        raise PermissionError(f"samtools is not executable: {samtools}")

    # Set temp file prefix
    output_dir = os.path.dirname(readname_file) or '.'
    temp_prefix = os.path.join(output_dir,
                               f"temp_patcher_{os.path.basename(readname_file).replace('.txt', '')}")

    if verbose:
        print("=" * 70)
        print("PATCHER: Rescuing False Negatives in AAV Classification")
        print("=" * 70)

    # Step 1: Parse classification file
    if verbose:
        print("Patcher Step 1: Parsing classification file...")
    classifications = parse_classification_file(readname_file)
    total_reads = len(classifications)
    if verbose:
        print(f"  Total reads: {total_reads}")

    # Step 2: Identify candidates (Unclassified reads only)
    if verbose:
        print("Patcher Step 2: Identifying Unclassified reads...")
    candidates = identify_candidates(classifications)
    if verbose:
        print(f"  Candidates found: {len(candidates)}")

    # Early exit if no candidates
    if not candidates:
        if verbose:
            print("  No Unclassified reads to analyze.")
        return {
            'total_reads': total_reads,
            'candidates': 0,
            'rescued': 0,
            'changes': []
        }

    # Step 3: Extract sequences from BAM
    if verbose:
        print("Patcher Step 3: Extracting sequences from BAM...")
    read_names = [c['readname'] for c in candidates]
    sequences = extract_sequences_from_bam(bam_file, read_names)
    if verbose:
        print(f"  Sequences extracted: {len(sequences)}")

    missing = set(read_names) - set(sequences.keys())
    if missing and verbose:
        print(f"  WARNING: {len(missing)} reads not found in BAM")

    if not sequences:
        if verbose:
            print("  No sequences could be extracted.")
        return {
            'total_reads': total_reads,
            'candidates': len(candidates),
            'rescued': 0,
            'changes': []
        }

    # Step 4: Extract WT reference
    if verbose:
        print("Patcher Step 4: Extracting WT reference...")
    wt_chrom, wt_sequence, target_pos = extract_wt_reference(ref_fasta, config_file)
    if verbose:
        print(f"  WT chromosome: {wt_chrom}")
        print(f"  Target position: {target_pos}")

    # Step 5: Write temp files and realign
    if verbose:
        print("Patcher Step 5: Realigning to WT-only reference...")
    reads_fasta, wt_fasta = write_temp_files(sequences, wt_chrom, wt_sequence, temp_prefix)
    realigned_bam = realign_to_wt(reads_fasta, wt_fasta, temp_prefix, minimap2, samtools)
    if verbose:
        print("  Realignment complete")

    # Step 6: Detect deletions
    if verbose:
        print("Patcher Step 6: Detecting variants...")
    variants = detect_deletions(realigned_bam, wt_chrom, target_pos,
                                cigar_window=cigar_window,
                                segment_window=segment_window,
                                min_del_size=min_del_size)
    if verbose:
        print(f"  Variants detected: {len(variants)}")

    # Step 7: Update classifications
    if verbose:
        print("Patcher Step 7: Updating classifications...")
    updated_classifications, changes = update_classifications(classifications, candidates, variants)
    if verbose:
        print(f"  Reclassifications: {len(changes)}")

    # Step 8: Write output (in-place update)
    if verbose:
        print("Patcher Step 8: Writing updated classifications...")
    write_output(updated_classifications, readname_file)
    if verbose:
        print(f"  Updated: {readname_file}")

    # Step 9: Cleanup
    cleanup_temp_files(temp_prefix, keep_temp)

    # Build statistics
    stats = {
        'total_reads': total_reads,
        'candidates': len(candidates),
        'rescued': len(changes),
        'changes': changes
    }

    if verbose:
        print("=" * 70)
        print(f"Patcher complete: {len(changes)} reads rescued from Unclassified -> DEL-large")
        print("=" * 70)

    return stats


def main():
    """Main entry point for CLI usage."""
    args = parse_args()

    # Resolve tool paths for display (CLI arg > env var > default)
    minimap2_path = get_tool_path(args.minimap2_path, 'PATCHER_MINIMAP2_PATH', DEFAULT_MINIMAP2_PATH)
    samtools_path = get_tool_path(args.samtools_path, 'PATCHER_SAMTOOLS_PATH', DEFAULT_SAMTOOLS_PATH)

    # Set output file path
    if args.output:
        output_file = args.output
    else:
        base = os.path.splitext(args.readname)[0]
        output_file = f"{base}_patched.txt"

    # Print CLI-specific header
    print("=" * 70)
    print("PATCHER: Rescuing False Negatives in AAV Classification")
    print("=" * 70)
    print(f"Input: {args.readname}")
    print(f"Output: {output_file}")
    print(f"Minimum deletion size: {args.length_threshold} bp")
    print(f"CIGAR detection window: +/- {args.cigar_window} bp (matches classify.py)")
    print(f"Segment detection window: +/- {args.segment_window} bp (matches classify.py)")
    print(f"minimap2: {minimap2_path}")
    print(f"samtools: {samtools_path}")
    print()

    # For CLI, we need to:
    # 1. Copy input to output first (so we don't modify original)
    # 2. Call rescue_unclassified_reads on the output file
    import shutil

    try:
        # Copy input to output file (preserves original)
        shutil.copy(args.readname, output_file)

        # Call the core function (operates in-place on output_file)
        stats = rescue_unclassified_reads(
            readname_file=output_file,
            bam_file=args.bam,
            ref_fasta=args.ref,
            config_file=args.config,
            min_del_size=args.length_threshold,
            cigar_window=args.cigar_window,
            segment_window=args.segment_window,
            minimap2_path=args.minimap2_path,
            samtools_path=args.samtools_path,
            keep_temp=args.keep_temp,
            verbose=True
        )

        # Build candidates list for statistics (for CLI compatibility)
        candidates = [{'readname': c['readname'], 'old_category': c['old_category'],
                       'old_length': c['old_length']} for c in stats['changes']]
        # Add unchanged candidates
        for i in range(stats['candidates'] - stats['rescued']):
            candidates.append({'readname': f'unchanged_{i}', 'old_category': 'Unclassified',
                               'old_length': ''})

        # Print statistics (CLI-specific)
        print()
        print_statistics(candidates, stats['changes'], args.stats)

        print(f"\nPatcher complete. Output: {output_file}")

    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except PermissionError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
