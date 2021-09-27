#!/usr/bin/env bash

in=$1   # FASTA to polish
r1=$2   # short reads (first in pairs)
r2=$3   # short reads (second in pairs)
rl=$4   # long reads
out=$5  # output prefix

# Make sure the necessary files exist.
if [[ ! -f "$in" ]] ; then echo $in" does not exist"; exit 1; fi
if [[ ! -f "$r1" ]] ; then echo $r1" does not exist"; exit 1; fi
if [[ ! -f "$r2" ]] ; then echo $r2" does not exist"; exit 1; fi
if [[ ! -f "$rl" ]] ; then echo $rl" does not exist"; exit 1; fi

printf "\n"
echo "PILON POLISHING WITH SHORT AND LONG READS"
echo "-----------------------------------------"
echo "Input FASTA: "$in
echo "Input short reads: "$r1", "$r2
echo "Input long reads: "$rl
printf "\n"

# Prepare the short-read alignments.
short_bam=pilon_short_"$RANDOM".bam
bwa index "$in"
bwa mem -t 16 "$in" "$r1" "$r2" | samtools sort > "$short_bam"
samtools index "$short_bam"

# Prepare the long-read alignments.
long_bam=pilon_long_"$RANDOM".bam
minimap2 -a -x map-ont -t 16 --secondary=no "$in" "$rl" | samtools sort > "$long_bam"
samtools index "$long_bam"

# Polish!
pilon --genome "$in" --frags "$short_bam" --nanopore "$long_bam" --output "$out"

# Make uppercase and one line per sequence.
temp=temp_"$RANDOM".fasta
seqtk seq "$out".fasta | awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' > "$temp"
mv "$temp" "$out".fasta

# Clean up.
rm "$short_bam" "$short_bam".bai
rm "$long_bam" "$long_bam".bai
rm "$in".amb "$in".ann "$in".bwt "$in".pac "$in".sa

printf "\n"
