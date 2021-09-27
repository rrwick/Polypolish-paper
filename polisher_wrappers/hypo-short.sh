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
echo "HYPO POLISHING WITH SHORT READS"
echo "-------------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# Prepare the short-read alignments.
bam=hypo_"$RANDOM".bam
bwa index "$in"
bwa mem -t 16 "$in" "$r1" "$r2" | samtools sort > "$bam"
samtools index "$bam"

# Get genome size and depth values.
genome_size=$(fast_count "$in" | cut -f3)
depth=$(samtools depth "$bam" | cut -f3 | awk '{total += $1} END {printf("%.0f\n", total/NR)}')

# Polish!
echo -e "$r1\n$r2" > read_filenames.txt
hypo -d "$in" -r @read_filenames.txt -s "$genome_size" -c "$depth" -b "$bam" -t 16 -o "$out".fasta

# Make uppercase and one line per sequence.
temp=temp_"$RANDOM".fasta
seqtk seq "$out".fasta | awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' > "$temp"
mv "$temp" "$out".fasta

# Clean up.
rm "$bam" "$bam".bai
rm read_filenames.txt
rm "$in".amb "$in".ann "$in".bwt "$in".pac "$in".sa
rm -r aux

printf "\n"



# This is the minimap2-based version of HyPo short-read polishing that I tried but ultimately didn't use:

# # Prepare the short-read alignments.
# bam=hypo_"$RANDOM".bam
# minimap2 --secondary=no --MD -ax sr -t 16 "$in" "$r1" "$r2" | samtools sort > "$bam"
# samtools index "$bam"

# # Get genome size and depth values.
# genome_size=$(fast_count "$in" | cut -f3)
# depth=$(samtools depth "$bam" | cut -f3 | awk '{total += $1} END {printf("%.0f\n", total/NR)}')

# # Polish!
# echo -e "$r1\n$r2" > read_filenames.txt
# hypo -d "$in" -r @read_filenames.txt -s "$genome_size" -c "$depth" -b "$bam" -t 16 -o "$out".fasta

# # Clean up.
# rm "$bam" "$bam".bai
# rm read_filenames.txt