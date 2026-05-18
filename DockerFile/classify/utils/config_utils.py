#!/usr/bin/env python3
"""
Configuration utilities for AAV Integration Classification Pipeline.

This module contains shared configuration parsing functions used by
classify_read.py and patcher.py to avoid circular imports.
"""


def read_aav_config(config_path):
    """
    Read in AAV vector config file and return a directory of AAV info.

    Config file should contain columns:
    AAV_Vector, Ref_Name, Ref_Type, ITR1_start, ITR1_end, gene_start, gene_end, ITR2_start, ITR2_end

    Args:
        config_path: Path to the AAV vector configuration TSV file

    Returns:
        dict: Nested dictionary with structure:
            {AAV_Vector: {Ref_Name: {field: value, ...}, ...}, ...}
            - Ref_Type is stored as string
            - Numeric fields are stored as int (or None if 'NA')
    """
    allvector = open(config_path, 'r').read().split('\n')
    allvector = [c for c in allvector if c != '']
    header = allvector[0][1:].split('\t')
    allvector = [c for c in allvector if c[0] != "#"]  # Remove header lines if there is any

    aavinfo = {}
    for line in allvector:
        line = line.split('\t')
        if line[0] not in aavinfo:
            aavinfo[line[0]] = {}
        aavinfo[line[0]][line[1]] = {}
        for i in range(len(header) - 2):
            if header[i + 2] == 'Ref_Type':
                aavinfo[line[0]][line[1]][header[i + 2]] = line[i + 2]
            else:
                if line[i + 2] != 'NA':
                    aavinfo[line[0]][line[1]][header[i + 2]] = int(line[i + 2])
                else:
                    aavinfo[line[0]][line[1]][header[i + 2]] = None
    return aavinfo
