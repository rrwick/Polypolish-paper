#!/usr/bin/env python3
"""
This script does pairwise global alignment of two FASTA files. The two files are assumed to have
one sequence each and to have already been normalised for strand and starting position.

It then prints the alignment on three lines:
* the first sequence
* the second sequence
* any differences (indicated with * characters)


Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import edlib
import gzip
import re
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit('Error: this script must be given two FASTA files as arguments')

    assembly_1_filename, assembly_2_filename = sys.argv[1], sys.argv[2]
    assembly_1, assembly_2 = load_fasta(assembly_1_filename), load_fasta(assembly_2_filename)

    if len(assembly_1) == 0 or len(assembly_2) == 0:
        print()
        sys.exit()

    if len(assembly_1) > 1:
        sys.exit(f'Error: {assembly_1_filename} has more than one sequence')
    if len(assembly_2) > 1:
        sys.exit(f'Error: {assembly_2_filename} has more than one sequence')

    seq_1, seq_2 = assembly_1[0], assembly_2[0]

    result = edlib.align(seq_1, seq_2, mode='NW', task='path')
    cigar = result['cigar']

    print_whole_alignment(assembly_1_filename, assembly_2_filename, seq_1, seq_2, cigar)


def print_whole_alignment(assembly_1_filename, assembly_2_filename, seq_1, seq_2, cigar):
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    aligned_seq_1, aligned_seq_2, differences = [], [], []
    for c in expanded_cigar:
        if c == '=':
            b_1 = seq_1[i]
            b_2 = seq_2[j]
            diff = ' '
            i += 1
            j += 1
            assert b_1 == b_2
        elif c == 'X':
            b_1 = seq_1[i]
            b_2 = seq_2[j]
            diff = '*'
            i += 1
            j += 1
            assert b_1 != b_2
        elif c == 'I':
            b_1 = seq_1[i]
            b_2 = '-'
            diff = '*'
            i += 1
        elif c == 'D':
            b_1 = '-'
            b_2 = seq_2[j]
            diff = '*'
            j += 1
        else:
            assert False
        aligned_seq_1.append(b_1)
        aligned_seq_2.append(b_2)
        differences.append(diff)
    assert i == len(seq_1)
    assert j == len(seq_2)
    aligned_seq_1 = ''.join(aligned_seq_1)
    aligned_seq_2 = ''.join(aligned_seq_2)
    differences = ''.join(differences)
    assert aligned_seq_1.replace('-', '') == seq_1
    assert aligned_seq_2.replace('-', '') == seq_2

    longer_name = max(len(assembly_1_filename), len(assembly_2_filename))
    assembly_1_filename = assembly_1_filename.rjust(longer_name)
    assembly_2_filename = assembly_2_filename.rjust(longer_name)

    print(f'{assembly_1_filename}: {aligned_seq_1}')
    print(f'{assembly_2_filename}: {aligned_seq_2}')
    spacer = ' ' * longer_name
    print(f'{spacer}  {differences}')


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_fasta(filename):
    try:
        fasta_seqs = []
        with get_open_func(filename)(filename, 'rt') as fasta_file:
            name = ''
            sequence = ''
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line[0] == '>':  # Header line = start of new contig
                    if name:
                        fasta_seqs.append(sequence)
                        sequence = ''
                    name = line[1:]
                else:
                    sequence += line
            if name:
                fasta_seqs.append(sequence)
        return fasta_seqs
    except FileNotFoundError:
        return []


if __name__ == '__main__':
    main()
