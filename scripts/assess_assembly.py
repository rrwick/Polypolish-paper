#!/usr/bin/env python3
"""
This script assesses the accuracy of an assembly in a number of ways:
 * alignment to a reference sequence
 * short-read-based assessment using BWA MEM and ALE
 * protein-based assessment with Prodigal

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
import gzip
import math
import pathlib
import re
import subprocess
import sys
import tempfile


def get_arguments():
    parser = argparse.ArgumentParser(description='Assembly assessor', add_help=False)

    required_args = parser.add_argument_group('Main arguments')
    required_args.add_argument('--header', action='store_true',
                               help='Print the header line and quit')
    required_args.add_argument('-g', '--genome', type=str,
                               help='The genome name being tested')
    required_args.add_argument('-n', '--name', type=str,
                               help='The name of the test')
    required_args.add_argument('-a', '--assembly', type=str,
                               help='Assembly to assess')
    required_args.add_argument('-r', '--reference', type=str,
                               help='Reference genome')
    required_args.add_argument('--repeats', type=str,
                               help='File defining repeat regions of the reference')

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
           'ref_filename', 'assembly_filename',
           'ref_length', 'assembly_length',
           'overall_errors', 'overall_identity', 'overall_qscore',
           'substitutions', 'insertions', 'deletions',
           'repeat_errors', 'repeat_identity', 'repeat_qscore',
           'non_repeat_errors', 'non_repeat_identity', 'non_repeat_qscore',
           'repeat_percentage', 'repeat_error_balance',
           'homopolymer_errors', 'non_homopolymer_errors',
           'lowest_100_bp_identity', 'lowest_1000_bp_identity',
           'prodigal_count', 'total_prodigal_length', 'mean_prodigal_length',
           'total_mapq', 'total_alignment_score',
           'ale_score', 'ale_place_avg', 'ale_insert_avg', 'ale_kmer_avg', 'ale_depth_score_avg']


def main():
    args = get_arguments()
    if args.header:
        print('\t'.join(HEADERS))
        sys.exit()

    if args.genome is None:
        sys.exit('Error: a value for --genome is required')
    if args.name is None:
        sys.exit('Error: a value for --name is required')

    assembly_seq, ref_seq = load_sequences(args)
    assembly_length, ref_length = len(assembly_seq), len(ref_seq)
    repeat_regions, total_repeat_size = load_repeats(args.repeats)

    expanded_cigar = align_assembly_to_ref(assembly_seq, ref_seq)

    errors, aligned_length, aligned_repeat_length, aligned_nonrepeat_length = \
        characterise_errors(expanded_cigar, assembly_seq, ref_seq, repeat_regions)

    overall_errors = count_errors(errors)
    repeat_errors = count_errors([e for e in errors if 'r' in e])
    nonrepeat_errors = count_errors([e for e in errors if 'r' not in e])
    assert repeat_errors + nonrepeat_errors == overall_errors

    substitutions = expanded_cigar.count('X')
    insertions = expanded_cigar.count('I')
    deletions = expanded_cigar.count('D')
    assert substitutions + insertions + deletions == overall_errors

    homopolymer_errors = count_errors([e for e in errors if 'h' in e])
    nonhomopolymer_errors = count_errors([e for e in errors if 'h' not in e])
    assert homopolymer_errors + nonhomopolymer_errors == overall_errors

    overall_identity = error_count_to_identity(overall_errors, aligned_length)
    repeat_identity = error_count_to_identity(repeat_errors, aligned_repeat_length)
    nonrepeat_identity = error_count_to_identity(nonrepeat_errors, aligned_nonrepeat_length)

    overall_qscore = identity_to_qscore(overall_identity)
    repeat_qscore = identity_to_qscore(repeat_identity)
    nonrepeat_qscore = identity_to_qscore(nonrepeat_identity)

    repeat_percentage = 100.0 * total_repeat_size / len(ref_seq)
    repeat_error_balance = get_repeat_error_balance(repeat_identity, nonrepeat_identity)

    lowest_100_bp_identity = get_lowest_window_identity(expanded_cigar, 100)
    lowest_1000_bp_identity = get_lowest_window_identity(expanded_cigar, 1000)

    prodigal_count, total_prodigal_length, mean_prodigal_length = \
        get_prodigal_scores(assembly_seq)

    total_mapq, total_alignment_score, \
        ale_score, place_score, insert_score, kmer_score, depth_score = \
        get_ale_scores(assembly_seq, args.short1, args.short2, args.threads)

    results = [args.genome, args.name,
               args.reference, args.assembly,
               str(ref_length), str(assembly_length),
               str(overall_errors), format_identity(overall_identity), overall_qscore,
               str(substitutions), str(insertions), str(deletions),
               str(repeat_errors), format_identity(repeat_identity), repeat_qscore,
               str(nonrepeat_errors), format_identity(nonrepeat_identity), nonrepeat_qscore,
               f'{repeat_percentage:.2f}%', repeat_error_balance,
               str(homopolymer_errors), str(nonhomopolymer_errors),
               f'{lowest_100_bp_identity:.0f}%', f'{lowest_1000_bp_identity:.1f}%',
               prodigal_count, total_prodigal_length, mean_prodigal_length,
               total_mapq, total_alignment_score,
               ale_score, place_score, insert_score, kmer_score, depth_score]

    print('\t'.join(results))
    print('', file=sys.stderr)


def load_sequences(args):
    print('', file=sys.stderr)
    if args.assembly is None:
        sys.exit('Error: an assembly filename is required')
    if args.reference is None:
        sys.exit('Error: a reference filename is required')
    
    print('Loading assembly sequence:  ', file=sys.stderr, end='')
    assembly = load_fasta(args.assembly)
    if len(assembly) == 0:
        sys.exit(f'\nError: no sequences in {args.assembly}')
    if len(assembly) > 1:
        sys.exit(f'\nError: {args.assembly} has more than one sequence')
    assembly_seq = assembly[0]
    print(f'{len(assembly_seq):,} bp', file=sys.stderr)

    print('Loading reference sequence: ', file=sys.stderr, end='')
    ref = load_fasta(args.reference)
    if len(ref) == 0:
        sys.exit(f'\nError: no sequences in {args.reference}')
    if len(ref) > 1:
        sys.exit(f'\nError: {args.reference} has more than one sequence')
    ref_seq = ref[0]
    print(f'{len(ref_seq):,} bp', file=sys.stderr)

    return assembly_seq, ref_seq


def load_repeats(repeat_filename):
    if repeat_filename is None:
        sys.exit('Error: no reference repeats supplied', file=sys.stderr)
    print('Loading reference repeats:  ', file=sys.stderr, end='')
    repeats = []
    total = 0
    with open(repeat_filename, 'rt') as repeat_file:
        for line in repeat_file:
            parts = line.strip().split('\t')
            if parts[0] == 'total':
                assert int(parts[1]) == total
            else:
                start, end = int(parts[0]), int(parts[1])
                repeats.append((start, end))
                total += (end - start)
    print(f'{len(repeats)} repeats loaded ({total} bp in total)', file=sys.stderr)
    return repeats, total


def align_assembly_to_ref(assembly_seq, ref_seq):
    print('', file=sys.stderr)
    print('Aligning assembly to reference:', file=sys.stderr)
    result = edlib.align(assembly_seq, ref_seq, mode='NW', task='path')
    cigar = result['cigar']
    expanded_cigar = get_expanded_cigar(cigar)
    
    matches = expanded_cigar.count('=')
    mismatches = expanded_cigar.count('X')
    insertions = expanded_cigar.count('I')
    deletions = expanded_cigar.count('D')
    
    print(f'  matches:    {matches:,}', file=sys.stderr)
    print(f'  mismatches: {mismatches:,}', file=sys.stderr)
    print(f'  insertions: {insertions:,}', file=sys.stderr)
    print(f'  deletions:  {deletions:,}', file=sys.stderr)
    print(f'  identity:   {100.0*matches/len(expanded_cigar):.5f}%', file=sys.stderr)

    return expanded_cigar


def characterise_errors(expanded_cigar, assembly_seq, ref_seq, repeat_regions):
    """
    This function returns a list of strings which characterise the errors in the assembly. They
    are in this format {ref_pos}-{cigar}, and they have an 'r' appended if they are in a repeat
    and an 'h' appended if they are a homopolymer indel error.
    """
    aligned_assembly_seq, aligned_ref_seq, aligned_repeat_regions, expanded_cigar = \
        align_seqs(expanded_cigar, assembly_seq, ref_seq, repeat_regions)

    errors = []
    e_pos, e_pos_ref, e_assembly, e_ref, e_cigar, e_repeat = None, None, None, None, None, None

    ref_i = 0
    for i, c in enumerate(expanded_cigar):
        if c == '=' and e_assembly is None:  # not in an error - nothing to do
            pass
        elif c == '=' and e_assembly is not None:  # finishing an error
            errors.append((e_pos, e_pos_ref, e_assembly, e_ref, e_cigar, e_repeat))
            e_pos, e_pos_ref, e_assembly, e_ref, e_cigar, e_repeat = \
                None, None, None, None, None, None
        elif c != '=' and e_assembly is None:  # starting an error
            e_pos = i
            e_pos_ref = ref_i
            e_assembly = aligned_assembly_seq[i]
            e_ref = aligned_ref_seq[i]
            e_cigar = c
            e_repeat = aligned_repeat_regions[i]
        elif c != '=' and e_assembly is not None:  # continuing an error
            e_assembly += aligned_assembly_seq[i]
            e_ref += aligned_ref_seq[i]
            e_cigar += c
            e_repeat += aligned_repeat_regions[i]
        if c == '=' or c == 'X' or c == 'D':
            ref_i += 1
    if e_assembly is not None:
        errors.append((e_pos, e_pos_ref, e_assembly, e_ref, e_cigar, e_repeat))
    assert ref_i == len(ref_seq)

    error_strings = []
    for e_pos, e_pos_ref, e_assembly, e_ref, e_cigar, e_repeat in errors:
        cigar_types = set(e_cigar)
        if cigar_types == {'D'} or cigar_types == {'I'}:  # simple indel
            homopolymer = is_indel_homopolymer(e_pos, len(e_cigar), aligned_assembly_seq,
                                               aligned_ref_seq)
        else:
            homopolymer = False
        error_string = f'{e_pos_ref}-{e_cigar}'
        if 'r' in e_repeat:
            error_string += 'r'
        if homopolymer:
            error_string += 'h'
        error_strings.append(error_string)

    aligned_length = len(aligned_repeat_regions)
    aligned_repeat_length = aligned_repeat_regions.count('r')
    aligned_nonrepeat_length = aligned_repeat_regions.count(' ')
    assert aligned_repeat_length + aligned_nonrepeat_length == aligned_length
    return error_strings, aligned_length, aligned_repeat_length, aligned_nonrepeat_length


def count_errors(error_strings):
    count = 0
    for e in error_strings:
        cigar = e.split('-')[1]
        cigar = cigar.replace('r', '').replace('h', '')
        assert set(cigar).issubset({'I', 'X', 'D'})
        count += len(cigar)
    return count


def align_seqs(expanded_cigar, assembly_seq, ref_seq, repeat_regions):
    aligned_assembly_seq, aligned_ref_seq, aligned_repeat_regions = [], [], []
    i, j = 0, 0
    for c in expanded_cigar:
        if is_in_repeat(j, repeat_regions):
            aligned_repeat_regions.append('r')
        else:
            aligned_repeat_regions.append(' ')
        if c == '=':
            assert assembly_seq[i] == ref_seq[j]
            aligned_assembly_seq.append(assembly_seq[i])
            aligned_ref_seq.append(ref_seq[j])
            i += 1
            j += 1
        elif c == 'X':
            assert assembly_seq[i] != ref_seq[j]
            aligned_assembly_seq.append(assembly_seq[i])
            aligned_ref_seq.append(ref_seq[j])
            i += 1
            j += 1
        elif c == 'I':
            aligned_assembly_seq.append(assembly_seq[i])
            aligned_ref_seq.append('-')
            i += 1
        elif c == 'D':
            aligned_assembly_seq.append('-')
            aligned_ref_seq.append(ref_seq[j])
            j += 1
    
    aligned_ref_seq = ''.join(aligned_ref_seq)
    aligned_assembly_seq = ''.join(aligned_assembly_seq)
    aligned_repeat_regions = ''.join(aligned_repeat_regions)
    assert len(aligned_ref_seq) == len(aligned_assembly_seq) == len(aligned_repeat_regions)

    return aligned_assembly_seq, aligned_ref_seq, aligned_repeat_regions, expanded_cigar


def error_count_to_identity(error_count, length):
    try:
        return 100.0 * (length - error_count) / length
    except ZeroDivisionError:
        return None


def format_identity(identity):
    if identity is None:
        return ''
    else:
        return f'{identity:.9f}%'


def identity_to_qscore(percent_identity):
    if percent_identity is None:
        return ''
    if percent_identity == 100.0:
        return 'inf'
    else:
        qscore = -10.0 * math.log10((100.0 - percent_identity) / 100.0)
        return f'{qscore:.2f}'


def get_repeat_error_balance(repeat_identity, nonrepeat_identity):
    if repeat_identity == 100.0 and nonrepeat_identity == 100.0:
        return ''
    if repeat_identity is None or nonrepeat_identity is None:
        return ''
    repeat_error_rate = 100.0 - repeat_identity
    nonrepeat_error_rate = 100.0 - nonrepeat_identity
    repeat_error_balance = repeat_error_rate / (repeat_error_rate + nonrepeat_error_rate)
    return f'{repeat_error_balance:.5f}'


def is_in_repeat(i, repeat_regions):
    for start, end in repeat_regions:
        if start <= i < end:
            return True
    return False


def get_lowest_window_identity(expanded_cigar, window_size):
    lowest_window_identity = 100.0
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


def is_indel_homopolymer(i, indel_size, seq_a, seq_b):
    """
    Checks to see if the indel in question is a homopolymer. In this function, the deletion can be
    on either of the two sequences.
    """
    if seq_a[i] == '-':
        for j in range(i, i+indel_size):
            assert seq_a[i] == '-'
            assert seq_b[i] != '-'
            return is_indel_homopolymer_2(i, indel_size, seq_a, seq_b)
    elif seq_b[i] == '-':
        for j in range(i, i+indel_size):
            assert seq_a[i] != '-'
            assert seq_b[i] == '-'
            return is_indel_homopolymer_2(i, indel_size, seq_b, seq_a)
    else:
        assert False


def is_indel_homopolymer_2(i, indel_size, seq_a, seq_b):
    """
    Checks to see if the indel in question is a homopolymer. In this function, the deletion is on
    sequence A.
    """
    bases = []
    for j in range(i, i+indel_size):
        assert seq_a[j] == '-'
        assert seq_b[j] != '-'
        bases.append(seq_b[j])

    # If there are multiple different bases in this indel, it can't be a homopolymer.
    if len(set(bases)) > 1:
        return False

    base = bases[0]
    assert all(b == base for b in bases)

    try:
        if seq_b[i-1] == base and seq_b[i-2] == base and seq_b[i-3] == base:
            return True
        j = i + indel_size - 1
        if seq_b[j+1] == base and seq_b[j+2] == base and seq_b[j+3] == base:
            return True
    except IndexError:
        pass

    return False


def get_prodigal_scores(assembly_seq):
    print('', file=sys.stderr)
    print('Getting proteins from Prodigal:', file=sys.stderr)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)

        assembly_filename = str(temp_dir / 'assembly.fasta')
        with open(assembly_filename, 'wt') as a:
            a.write(f'>assembly\n{assembly_seq}\n')

        prodigal_out = str(temp_dir / 'proteins.fasta')
        prodigal_command = f'prodigal -a {prodigal_out} -q -i {assembly_filename}'
        print('  ' + prodigal_command, file=sys.stderr)
        subprocess.run(prodigal_command, shell=True, capture_output=True)

        proteins = load_fasta(prodigal_out)
        prodigal_count = len(proteins)
        total_prodigal_length = sum(len(p) for p in proteins)
        mean_prodigal_length = total_prodigal_length / prodigal_count

        prodigal_count = str(prodigal_count)
        total_prodigal_length = str(total_prodigal_length)
        mean_prodigal_length = f'{mean_prodigal_length:.3f}'

        print(f'  Protein count:        {prodigal_count}', file=sys.stderr)
        print(f'  Total protein length: {total_prodigal_length}', file=sys.stderr)
        print(f'  Mean protein length:  {mean_prodigal_length}', file=sys.stderr)

    return prodigal_count, total_prodigal_length, mean_prodigal_length


def get_ale_scores(assembly_seq, reads_1, reads_2, threads):
    if reads_1 is None or reads_2 is None:
        return '', '', '', '', '', '', ''

    print('', file=sys.stderr)
    print('Getting ALE scores:', file=sys.stderr)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)

        assembly_filename = str(temp_dir / 'assembly.fasta')
        with open(assembly_filename, 'wt') as a:
            a.write(f'>assembly\n{assembly_seq}\n')

        sam = str(temp_dir / 'alignments.sam')
        ale_out = str(temp_dir / 'ale.txt')

        index_command = f'bwa index {assembly_filename}'
        print('  ' + index_command, file=sys.stderr)
        subprocess.run(index_command, shell=True, capture_output=True)

        align_command = f'bwa mem -t {threads} {assembly_filename} {reads_1} {reads_2} > {sam}'
        print('  ' + align_command, file=sys.stderr)
        subprocess.run(align_command, shell=True, capture_output=True)

        total_mapq, total_alignment_score = get_total_mapq_and_alignment_score(sam)
        print(f'  Total MAPQ:            {total_mapq}', file=sys.stderr)
        print(f'  Total alignment score: {total_alignment_score}', file=sys.stderr)
        total_mapq, total_alignment_score = str(total_mapq), str(total_alignment_score)

        ale_command = f'ALE {sam} {assembly_filename} {ale_out}'
        print('  ' + ale_command, file=sys.stderr)
        subprocess.run(ale_command, shell=True, capture_output=True)

        ale_score, place_score, insert_score, kmer_score, depth_score = None, None, None, None, None
        with open(ale_out, 'rt') as a:
            for line in a:
                line = line.strip()
                if line.startswith('# ALE_score: '):
                    ale_score = line.split('# ALE_score: ')[1]
                if line.startswith('# placeAvg: '):
                    place_score = line.split('# placeAvg: ')[1]
                if line.startswith('# insertAvg: '):
                    insert_score = line.split('# insertAvg: ')[1]
                if line.startswith('# kmerAvg: '):
                    kmer_score = line.split('# kmerAvg: ')[1]
                if line.startswith('# depthScoreAvg: '):
                    depth_score = line.split('# depthScoreAvg: ')[1]
        print(f'  ALE score:       {ale_score}', file=sys.stderr)
        print(f'  placement score: {place_score}', file=sys.stderr)
        print(f'  insert score:    {insert_score}', file=sys.stderr)
        print(f'  k-mer score:     {kmer_score}', file=sys.stderr)
        print(f'  depth score:     {depth_score}', file=sys.stderr)

    return total_mapq, total_alignment_score, \
        ale_score, place_score, insert_score, kmer_score, depth_score


def get_total_mapq_and_alignment_score(sam):
    total_mapq, total_alignment_score = 0, 0
    with open(sam, 'rt') as s:
        for line in s:
            if line.startswith('@'):  # header
                continue
            parts = line.strip().split('\t')
            mapq = int(parts[4])
            if mapq == 255:
                mapq = 0
            alignment_score = 0
            for p in parts:
                if p.startswith('AS:i:'):
                    alignment_score = int(p[5:])
            total_mapq += mapq
            total_alignment_score += alignment_score
    return total_mapq, total_alignment_score


if __name__ == '__main__':
    main()
