#!/usr/bin/env python3
"""
This script takes in a bunch of verbose time reports (from /usr/bin/time -v) and outputs the
mean/median time and memory usage.

Copyright 2021 Ryan Wick (rrwick@gmail.com)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import statistics
import sys


def main():
    filenames = sys.argv[1:]
    time_taken, memory_used = [], []
    for filename in filenames:
        time_taken.append(get_time_taken(filename))
        memory_used.append(get_memory_used(filename))
    
    time_taken = [t for t in time_taken if t is not None]
    memory_used = [m for m in memory_used if m is not None]
    assert len(time_taken) == len(memory_used)

    total_hours = sum(time_taken) / 3600.0
    mean_minutes = statistics.mean(time_taken) / 60.0
    mean_gigabytes = statistics.mean(memory_used) / 1000000
    median_minutes = statistics.median(time_taken) / 60.0
    median_gigabytes = statistics.median(memory_used) / 1000000
    
    print()
    print(f'From {len(time_taken)} files:')
    print()
    print(f'total time (hours):        {total_hours:.2f}')
    print()
    print(f'mean time (minutes):       {mean_minutes:.2f}')
    print(f'mean memory (gigabytes):   {mean_gigabytes:.2f}')
    print()
    print(f'median time (minutes):     {median_minutes:.2f}')
    print(f'median memory (gigabytes): {median_gigabytes:.2f}')
    print()


def get_time_taken(filename):
    with open(filename, 'rt') as f:
        for line in f:
            if 'Elapsed (wall clock) time (h:mm:ss or m:ss): ' in line:
                time_str = line.strip().split('(h:mm:ss or m:ss): ')[1]
                return convert_time_to_seconds(time_str)
    return None


def convert_time_to_seconds(time_str):
    parts = time_str.split(':')
    if len(parts) == 2:
        hours = 0.0
        minutes = float(parts[0])
        seconds = float(parts[1])
    elif len(parts) == 3:
        hours = float(parts[0])
        minutes = float(parts[1])
        seconds = float(parts[2])
    else:
        sys.exit(f'Error: could not parse time in {time_str}')
    return 3600.0 * hours + 60.0 * minutes + seconds


def get_memory_used(filename):
    with open(filename, 'rt') as f:
        for line in f:
            if 'Maximum resident set size (kbytes): ' in line:
                return int(line.strip().split('(kbytes): ')[1])
    return None


if __name__ == '__main__':
    main()
