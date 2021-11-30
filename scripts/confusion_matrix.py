#!/usr/bin/env python3
"""
This script reports the following counts for a polished assembly:
* True positives:  positions where the unpolished assembly is wrong and the polished assembly is
                   right (fixed error)
* True negatives:  positions where the unpolished assembly is right and the polished assembly is
                   still right (correctly left alone)
* False positives: positions where the unpolished assembly is right and the polished assembly is
                   wrong (introduced error)
* False negatives: positions where the unpolished assembly is wrong and the polished assembly is
                   wrong (unfixed error) - includes cases where one error is replaced by a
                   different error

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
import gzip
import pathlib
import subprocess
import sys
import tempfile


def get_arguments():
    parser = argparse.ArgumentParser(description='Polishing confusion matrix', add_help=False)

    required_args = parser.add_argument_group('Main arguments')
    required_args.add_argument('-r', '--reference', type=str,
                               help='Reference genome')
    required_args.add_argument('-u', '--unpolished', type=str,
                               help='Unpolished genome')
    required_args.add_argument('-p', '--polished', type=str,
                               help='Polished genome')
    required_args.add_argument('--header', action='store_true',
                               help='Print the header line and quit')
    required_args.add_argument('-g', '--genome', type=str,
                               help='The genome name being tested')
    required_args.add_argument('-n', '--name', type=str,
                               help='The name of the test')

    read_args = parser.add_argument_group('Reads')
    read_args.add_argument('-1', '--short1', type=str,
                           help='Illumina read file (first reads in pair)')
    read_args.add_argument('-2', '--short2', type=str,
                           help='Illumina read file (second reads in pair)')

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('--threads', type=int, default=16,
                            help='Number of CPU threads for alignment')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')

    args = parser.parse_args()
    return args


HEADERS = ['genome', 'test_name',
           'true_positives', 'true_negatives',
           'false_positives', 'false_negatives']


def main():
    args = get_arguments()
    if args.header:
        print('\t'.join(HEADERS))
        sys.exit()
    check_requirements(args)
    reference_seq, unpolished_seq, polished_seq, ref_length = \
        get_msa(args.reference, args.unpolished, args.polished)
    tp, tn, fp, fn = count_positions(reference_seq, unpolished_seq, polished_seq)
    assert tp + tn + fp + fn == ref_length
    print(f'{args.genome}\t{args.name}\t{tp}\t{tn}\t{fp}\t{fn}')


def check_requirements(args):
    if args.genome is None:
        sys.exit('Error: a value for --genome is required')
    if args.name is None:
        sys.exit('Error: a value for --name is required')
    if args.reference is None:
        sys.exit('Error: a value for --reference is required')
    if args.unpolished is None:
        sys.exit('Error: a value for --unpolished is required')
    if args.polished is None:
        sys.exit('Error: a value for --polished is required')


def get_msa(reference, unpolished, polished):
    reference_seq, unpolished_seq, polished_seq = load_sequences(reference, unpolished, polished)
    ref_length = len(reference_seq)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        all_seqs_filename = temp_dir / '2_all_seqs.fasta'
        msa_filename = temp_dir / '3_msa.fasta'
        with open(all_seqs_filename, 'wt') as all_seqs:
            all_seqs.write(f'>reference\n{reference_seq}\n')
            all_seqs.write(f'>unpolished\n{unpolished_seq}\n')
            all_seqs.write(f'>polished\n{polished_seq}\n')
        trycycler_command = ['trycycler', 'msa', '-c', str(temp_dir)]

        process = subprocess.Popen(trycycler_command, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, encoding='utf-8')
        while True:
            realtime_output = process.stdout.readline()
            if realtime_output == '' and process.poll() is not None:
                break
            if realtime_output and not realtime_output.startswith('pieces:'):
                print(realtime_output.strip(), file=sys.stderr, flush=True)

        msa_seqs = load_fasta_with_names(msa_filename)
        assert len(msa_seqs) == 3
        assert 'reference' in msa_seqs
        assert 'unpolished' in msa_seqs
        assert 'polished' in msa_seqs

        return msa_seqs['reference'], msa_seqs['unpolished'], msa_seqs['polished'], ref_length


def count_positions(reference_seq, unpolished_seq, polished_seq):
    tp, tn, fp, fn = 0, 0, 0, 0
    for start, end in iterate_seq(reference_seq):
        reference_slice = reference_seq[start:end]
        unpolished_slice = unpolished_seq[start:end]
        polished_slice = polished_seq[start:end]

        unpolished_correct = (reference_slice == unpolished_slice)
        polished_correct = (reference_slice == polished_slice)

        if unpolished_correct:
            if polished_correct:
                tn += 1
            else:
                fp += 1
        else:
            if polished_correct:
                tp += 1
            else:
                fn += 1

    return tp, tn, fp, fn


def iterate_seq(seq):
    """
    This generator iterates through a sequence one base at a time, returning start/end slice
    positions. Gaps are included after a base, so each returned slice positions will include
    exactly one base followed by zero or more gaps.
    """
    if not seq:
        return
    start = 0
    while True:
        end = start + 1
        while True:
            if end < len(seq) and seq[end] == '-':
                end += 1
            else:
                break
        yield start, end
        start = end
        if start >= len(seq):
            return



def load_sequences(reference_filename, unpolished_filename, polished_filename):
    print('Loading reference sequence:  ', file=sys.stderr, end='')
    reference = load_fasta(reference_filename)
    if len(reference) == 0:
        sys.exit(f'\nError: no sequences in {reference_filename}')
    if len(reference) > 1:
        sys.exit(f'\nError: {reference_filename} has more than one sequence')
    reference_seq = reference[0]
    print(f'{len(reference_seq):,} bp', file=sys.stderr)

    print('Loading unpolished sequence: ', file=sys.stderr, end='')
    unpolished = load_fasta(unpolished_filename)
    if len(unpolished) == 0:
        sys.exit(f'\nError: no sequences in {unpolished_filename}')
    if len(unpolished) > 1:
        sys.exit(f'\nError: {unpolished_filename} has more than one sequence')
    unpolished_seq = unpolished[0]
    print(f'{len(unpolished_seq):,} bp', file=sys.stderr)

    print('Loading polished sequence:   ', file=sys.stderr, end='')
    polished = load_fasta(polished_filename)
    if len(polished) == 0:
        sys.exit(f'\nError: no sequences in {polished_filename}')
    if len(polished) > 1:
        sys.exit(f'\nError: {polished_filename} has more than one sequence')
    polished_seq = polished[0]
    print(f'{len(polished_seq):,} bp', file=sys.stderr)

    return reference_seq, unpolished_seq, polished_seq


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
                if line[0] == '>':  # header line = start of new contig
                    if name:
                        fasta_seqs.append(sequence)
                        sequence = ''
                    name = line[1:]
                else:
                    sequence += line.upper()
            if name:
                fasta_seqs.append(sequence)
        return fasta_seqs
    except FileNotFoundError:
        return []


def load_fasta_with_names(filename):
    fasta_seqs = {}
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs[name.split()[0]] = ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            fasta_seqs[name.split()[0]] = ''.join(sequence)
    return fasta_seqs


if __name__ == '__main__':
    main()
