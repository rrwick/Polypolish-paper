#!/usr/bin/env bash

in=$1   # FASTA to polish
r1=$2   # short reads (first in pairs)
r2=$3   # short reads (second in pairs)
out=$4  # output prefix

# Make sure the necessary files exist.
if [[ ! -f "$in" ]] ; then echo $in" does not exist"; exit 1; fi
if [[ ! -f "$r1" ]] ; then echo $r1" does not exist"; exit 1; fi
if [[ ! -f "$r2" ]] ; then echo $r2" does not exist"; exit 1; fi

printf "\n"
echo "RACON POLISHING WITH SHORT READS"
echo "--------------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# Prepare a single file of short reads.
combined_reads=combined_reads_"$RANDOM".fastq
gunzip -c "$r1" > "$combined_reads"
gunzip -c "$r2" >> "$combined_reads"

# Prepare the short-read alignments.
sam=racon_"$RANDOM".sam
minimap2 -t 16 -a -x sr "$in" "$combined_reads" > "$sam"

# Polish!
racon -t 16 --no-trimming "$combined_reads" "$sam" "$in" > "$out".fasta

# Make uppercase and one line per sequence.
temp=temp_"$RANDOM".fasta
seqtk seq "$out".fasta | awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' > "$temp"
mv "$temp" "$out".fasta

# Clean up.
rm "$sam" "$combined_reads"

printf "\n"
