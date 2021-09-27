#!/usr/bin/env python3
"""
This script takes a genome as input and outputs (to stdout) a version of the genome with errors
added. The errors can be a few different types:
* small errors: single-bp substitutions and indels
* large errors: multi-bp substitutions and indels
* homopolymer errors: changes in homopolymer length

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
import collections
import gzip
import numpy
import random
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Add small errors to genome')

    parser.add_argument('input_genome', type=str,
                        help='FASTA file of genome to which errors will be added')
    parser.add_argument('error_type', type=str,
                        choices=['small', 'large', 'homopolymer'],
                        help='Type of errors to be added')
    parser.add_argument('--error_rate', type=float, default=0.0001,
                        help='Target error rate')

    large_args = parser.add_argument_group('Large error options',
                                           description='Options for when error_type = large')
    large_args.add_argument('--large_error_size', type=float, default=5.0,
                            help='Mean error size for large errors')
    large_args.add_argument('--max_large_error_size', type=int, default=25,
                            help='Max error size for large errors')

    homo_args = parser.add_argument_group('Homopolymer error options',
                                          description='Options for when error_type = homopolymer')
    homo_args.add_argument('--min_homopolymer_len', type=int, default=3,
                           help='Only homopolymers this length or longer will be changed')
    homo_args.add_argument('--insertions_only', action='store_true',
                           help='Only make homopolymers longer (not shorter)')

    args = parser.parse_args()
    return args


def main():
    random.seed(0)
    args = get_arguments()
    print(f'\nAdding {args.error_type} errors to {args.input_genome}', file=sys.stderr, flush=True)
    for name, seq in load_fasta(args.input_genome):
        if args.error_type == 'small':
            seq_with_errors = add_small_errors(seq, args.error_rate)
            print(f'>{name}_small_errors')
        elif args.error_type == 'large':
            seq_with_errors = add_large_errors(seq, args.error_rate, args.large_error_size,
                                               args.max_large_error_size)
            print(f'>{name}_large_errors')
        elif args.error_type == 'homopolymer':
            seq_with_errors = add_homopolymer_errors(seq, args.error_rate, args.insertions_only,
                                                     args.min_homopolymer_len)
            print(f'>{name}_homopolymer_errors')
        else:
            assert False
        print(seq_with_errors)
    print(f'Done!\n', file=sys.stderr)


def add_small_errors(seq, error_rate):
    seq_with_errors = [b for b in seq]
    positions = list(range(0, len(seq)))
    random.shuffle(positions)
    target_errors = int(round(error_rate * len(seq)))
    errors_added = 0
    for i in positions:
        base = seq_with_errors[i]
        error_type = random.randint(0, 2)
        if error_type == 0:  # substitution
            seq_with_errors[i] = get_random_different_base(base)
        elif error_type == 1:  # insertion
            if random.randint(0, 1) == 0:  # insertion before base
                seq_with_errors[i] = get_random_base() + base
            else:  # insertion after base
                seq_with_errors[i] = base + get_random_base()
        elif error_type == 2:  # deletion
            seq_with_errors[i] = ''
        else:
            assert False
        errors_added += 1
        if errors_added >= target_errors:
            break
    return ''.join(seq_with_errors)


def add_large_errors(seq, error_rate, mean_error_size, max_error_size):
    seq_with_errors = [b for b in seq]
    positions = list(range(0, len(seq)))
    random.shuffle(positions)
    target_errors = int(round(error_rate * len(seq)))
    errors_added = 0
    for i in positions:
        error_size = numpy.random.geometric(1.0 / mean_error_size)
        while error_size > max_error_size:
            error_size = numpy.random.geometric(1.0 / mean_error_size)
        error_type = random.randint(0, 2)
        if error_type == 0:  # substitution
            errors_added += int(round(error_size * 0.85))  # adjust down to compensate for alignment
            for _ in range(error_size):
                seq_with_errors[i] = get_random_different_base(seq_with_errors[i])
                i += 1
                if i == len(seq):
                    break
        elif error_type == 1:  # insertion
            errors_added += error_size
            new_seq = ''.join([get_random_base() for _ in range(error_size)])
            if random.randint(0, 1) == 0:  # insertion before base
                seq_with_errors[i] = new_seq + seq_with_errors[i]
            else:  # insertion after base
                seq_with_errors[i] = seq_with_errors[i] + new_seq
        elif error_type == 2:  # deletion
            errors_added += error_size
            for _ in range(error_size):
                seq_with_errors[i] = ''
                i += 1
                if i == len(seq):
                    break
        else:
            assert False
        if errors_added >= target_errors:
            break
    return ''.join(seq_with_errors)


def add_homopolymer_errors(seq, error_rate, insertions_only, min_len):
    target_errors = int(round(error_rate * len(seq)))
    target_min = target_errors * 0.99
    target_max = target_errors * 1.01
    stdev = 0.025
    partitioned_seq = partition_seq_by_homopolymers(seq)
    while True:
        errors_added = 0
        seq_with_errors = []
        for piece in partitioned_seq:
            if len(piece) < min_len:  # not a homopolymer
                seq_with_errors.append(piece)
            else:
                assert len(set(piece)) == 1
                base = piece[0]
                old_length = len(piece)
                while True:
                    new_length = int(round(numpy.random.normal(old_length, old_length * stdev)))
                    if new_length >= old_length or not insertions_only:
                        break
                if insertions_only:
                    assert new_length >= old_length
                if new_length < 1:
                    new_length = 1
                seq_with_errors.append(base * new_length)
                errors_added += abs(new_length - old_length)
        print(f'stdev = {stdev:.5f}, {errors_added} / {target_errors} errors added',
              file=sys.stderr)
        if errors_added < target_min:
            stdev += 0.0005
        elif errors_added > target_max:
            stdev -= 0.0005
        else:
            return ''.join(seq_with_errors)


def partition_seq_by_homopolymers(seq):
    """
    This function takes in a sequence and returns in a list form where homopolymers are grouped.
    For example:
      * input: 'CGAAAACACCGACGGGT'
      * output: ['C', 'G, 'AAAA, 'C, 'A, 'CC, 'G, 'A, 'C, 'GGG, 'T']
    """
    partitioned_seq = []
    current_piece = ''
    for base in seq:
        if not current_piece:  # first base
            current_piece = base
        elif base == current_piece[0]:  # extending a homopolymer
            current_piece += base
        else:  # not extending a homopolymer
            partitioned_seq.append(current_piece)
            current_piece = base
    partitioned_seq.append(current_piece)
    assert ''.join(partitioned_seq) == seq
    return partitioned_seq


def load_fasta(fasta_filename, include_full_header=False):
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
                    if include_full_header:
                        fasta_seqs.append((name.split()[0], name, ''.join(sequence)))
                    else:
                        fasta_seqs.append((name.split()[0], ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            if include_full_header:
                fasta_seqs.append((name.split()[0], name, ''.join(sequence)))
            else:
                fasta_seqs.append((name.split()[0], ''.join(sequence)))
    return fasta_seqs


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
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


def get_random_base():
    return {0: 'A', 1: 'C', 2: 'G', 3: 'T'}[random.randint(0, 3)]


def get_random_different_base(base):
    while True:
        new_base = get_random_base()
        if new_base != base:
            return new_base


if __name__ == '__main__':
    main()
