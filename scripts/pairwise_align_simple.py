#!/usr/bin/env python3
"""
This script does pairwise global alignment of two FASTA files. The two files are assumed to have
one sequence each and to have already been normalised for strand and starting position.

It then outputs:
* the overall sequence identity of the alignment and the corresponding qscore
* the lowest sequence identity of a 100 bp sliding window
* the lowest sequence identity of a 1000 bp sliding window
* the total number of single-bp errors in the alignment (e.g. a 5-bp indel add 5 to this total)


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
import math
import re
import sys


def main():
    if len(sys.argv) == 1:
        print_header()
        sys.exit(0)
    if len(sys.argv) != 3:
        sys.exit('Error: this script must be given two FASTA files as arguments')

    assembly_1_filename, assembly_2_filename = sys.argv[1], sys.argv[2]
    assembly_1, assembly_2 = load_fasta(assembly_1_filename), load_fasta(assembly_2_filename)

    if len(assembly_1) == 0 or len(assembly_2) == 0:
        print()
        sys.exit()

    assembly_1 = sorted(assembly_1, key=lambda x: len(x), reverse=True)
    assembly_2 = sorted(assembly_2, key=lambda x: len(x), reverse=True)

    seq_1, seq_2 = assembly_1[0], assembly_2[0]

    result = edlib.align(seq_1, seq_2, mode='NW', task='path')
    cigar = result['cigar']

    errors, identity = get_alignment_stats(cigar)
    qscore = identity_to_qscore(identity)
    lowest_100_bp_identity = get_lowest_window_identity(cigar, 100)
    lowest_1000_bp_identity = get_lowest_window_identity(cigar, 1000)

    print(f'{assembly_1_filename}\t{assembly_2_filename}\t'
          f'{errors}\t{identity:.9f}%\t{qscore}\t'
          f'{lowest_100_bp_identity:.1f}%\t{lowest_1000_bp_identity:.1f}%')


def print_header():
    print('assembly_1_filename\tassembly_2_filename\t'
          'errors\tidentity\tqscore\t'
          'lowest_100_bp_identity\tlowest_1000_bp_identity')


def get_alignment_stats(cigar):
    expanded_cigar = get_expanded_cigar(cigar)
    length = 0
    errors = 0
    for c in expanded_cigar:
        if c != '=':
            errors += 1
        length += 1
    identity = 100.0 * (length - errors) / length

    return errors, identity


def identity_to_qscore(percent_identity):
    if percent_identity == 100.0:
        return 'inf'
    else:
        qscore = -10.0 * math.log10((100.0 - percent_identity) / 100.0)
        return f'{qscore:.2f}'


def get_lowest_window_identity(cigar, window_size):
    lowest_window_identity = 100.0
    expanded_cigar = get_expanded_cigar(cigar)
    for i in range(0, len(expanded_cigar) - window_size):
        cigar_slice = expanded_cigar[i:i+window_size]
        window_identity = 100.0 * cigar_slice.count('=') / window_size
        if window_identity < lowest_window_identity:
            lowest_window_identity = window_identity
    return lowest_window_identity


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
