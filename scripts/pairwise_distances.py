#!/usr/bin/env python3

import argparse
import edlib
import re
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Find all pairwise distances')

    parser.add_argument('fasta_files', type=str, nargs='+',
                        help='FASTA file with all sequences to be compared')

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
