#!/usr/bin/env python3
"""
This script reads SAM alignments from stdin and output the regions of the genome which are
repetitive. Whenever one read aligns to two or more places in the genome, all of the genomic
regions under those alignments are considered to be repetitive. I.e. in this context, 'repetitive'
is defined with respect to the read alignments. So shorter reads will result in more repetitive
regions than longer reads.

Usage:
  bwa mem -t 16 -a genome.fasta reads.fastq | samtools view -F 2048 | find_repetitve_regions.py 

Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import collections
import fileinput
import re
import sys


def main():
    alignments_by_read = get_alignments()
    repeat_bases = get_repeat_bases(alignments_by_read)
    total_repeat_size = sum(len(r) for r in repeat_bases.values())
    repeat_ranges = get_repeat_ranges(repeat_bases)
    for reference_name, ranges in repeat_ranges.items():
        for start, end in ranges:
            print(f'{reference_name}\t{start}\t{end}')
    print(f'total\t{total_repeat_size}')


def get_alignments():
    """
    This function reads SAM alignments from stdin. It returns the alignments (with only their
    start position and CIGAR) grouped by read.
    """
    alignments_by_read = collections.defaultdict(list)
    for line in fileinput.input():
        if line.startswith('@'):  # header section
            continue
        line = line.strip()
        parts = line.split('\t')
        read_name = parts[0]
        reference_name = parts[2]
        start_pos = int(parts[3]) - 1
        cigar = parts[5]
        alignments_by_read[read_name].append((reference_name, start_pos, cigar))
    return alignments_by_read


def get_repeat_bases(alignments_by_read):
    """
    This function returns any bases in the reference which are repeats. It looks for any reads
    with multiple alignments, and considers the bases covered by these alignments to be repetitive.
    """
    repeat_bases = collections.defaultdict(set)
    for read_name, alignments in alignments_by_read.items():
        if len(alignments) < 2:
            continue
        for reference_name, start_pos, cigar in alignments:
            end_pos = get_end_pos(start_pos, cigar)
            for i in range(start_pos, end_pos):
                repeat_bases[reference_name].add(i)
    return repeat_bases


def cigar_starts_and_ends_with_match(cigar):
    cigar_parts = re.findall(r'\d+[MIDNSH=X]', cigar)
    first_part, last_part = cigar_parts[0], cigar_parts[-1]
    return first_part[-1] == 'M' and last_part[-1] == 'M'


def get_end_pos(start_pos, cigar):
    """
    Given a start position and a CIGAR from a SAM line, this function returns the end position
    of the alignment. It adds to the start position for each CIGAR operation which 'consumes' a
    reference base.
    """
    cigar_parts = re.findall(r'\d+[MIDNSH=X]', cigar)
    pos = start_pos
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        if letter == 'M' or letter == 'D' or letter == 'N' or letter == '=' or letter == 'X':
            pos += size
    return pos


def get_repeat_ranges(repeat_bases):
    ranges = {}
    for reference_name, bases in repeat_bases.items():
        ranges[reference_name] = get_repeat_ranges_one_ref(bases)
    return ranges


def get_repeat_ranges_one_ref(repeat_bases):
    """
    This function takes in the set of repeat bases and returns the discrete integer ranges they
    cover.
    """
    ranges = []
    current_range = None
    for b in sorted(repeat_bases):
        if current_range is None:    # first range
            current_range = [b, b+1]
        elif b == current_range[1]:  # extending a range
            current_range[1] += 1
        else:                        # end one range and start a new one
            ranges.append(current_range)
            current_range = [b, b+1]
    if current_range is not None:
        ranges.append(current_range)
    return ranges


if __name__ == '__main__':
    main()
