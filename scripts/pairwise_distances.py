#!/usr/bin/env python3
"""
This script takes a group of FASTA files as input. It loads the longest sequence from each (assumed
to be the chromosome) and performs all pairwise alignments between them, ouputting a PHYLIP matrix
and printing the total distance to stderr.

Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import argparse
import edlib
import re
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Find all pairwise distances')

    parser.add_argument('fasta_files', type=str, nargs='+',
                        help='FASTA files with sequences to be compared')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    seqs = [load_longest_seq(f) for f in args.fasta_files]
    seqs = [s for s in seqs if s[0] is not None]

    print(len(seqs))
    distances = {}
    total_distance = 0
    for i in range(len(seqs)):
        name_a, seq_a = seqs[i]
        a_distances = []
        for j in range(len(seqs)):
            name_b, seq_b = seqs[j]
            if (name_a, name_b) in distances:
                distance = distances[(name_a, name_b)]
            else:
                distance = get_distance(seq_a, seq_b)
                distances[(name_a, name_b)] = distance
                distances[(name_b, name_a)] = distance
                total_distance += distance
            a_distances.append(distance)
        a_distances = '\t'.join([str(d) for d in a_distances])
        print(f'{name_a}\t{a_distances}')
    print(total_distance, file=sys.stderr)


def load_longest_seq(fasta_filename):
    longest_seq = '', ''
    try:
        for _, seq in load_fasta(fasta_filename):
            if len(seq) > len(longest_seq):
                longest_seq = seq
    except FileNotFoundError:
        return None, None
    return fasta_filename, longest_seq


def get_distance(seq_a, seq_b):
    assert len(seq_a) > 1000000 and len(seq_b) > 1000000  # make sure we have the chromosome
    if seq_a == seq_b:
        return 0
    else:
        return get_edlib_all_distance(seq_a, seq_b)


def get_edlib_all_distance(seq_a, seq_b):
    result = edlib.align(seq_a, seq_b, mode="NW", task="path")
    cigar = result['cigar']
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    distance = 0
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        if letter != '=':
            distance += size
    return distance


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


def load_fasta(fasta_filename):
    if get_compression_type(fasta_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    fasta_seqs = []
    with open_func(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            fasta_seqs.append((name.split()[0], ''.join(sequence)))
    return fasta_seqs


if __name__ == '__main__':
    main()
