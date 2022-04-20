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
echo "POLYPOLISH POLISHING WITH SHORT READS"
echo "-------------------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# Prepare the short-read alignments.
sam1=polypolish_1_"$RANDOM".sam
sam2=polypolish_2_"$RANDOM".sam
filtered_sam1=polypolish_1_filtered_"$RANDOM".sam
filtered_sam2=polypolish_2_filtered_"$RANDOM".sam
bwa index "$in"
bwa mem -t 16 -a "$in" "$r1" > "$sam1"
bwa mem -t 16 -a "$in" "$r2" > "$sam2"
polypolish_insert_filter --in1 "$sam1" --in2 "$sam2" --out1 "$filtered_sam1" --out2 "$filtered_sam2"

# Polish!
polypolish "$in" "$filtered_sam1" "$filtered_sam2" > "$out".fasta

# Make uppercase and one line per sequence.
temp=temp_"$RANDOM".fasta
seqtk seq "$out".fasta | awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' > "$temp"
mv "$temp" "$out".fasta

# Clean up.
rm "$sam1" "$sam2" "$filtered_sam1" "$filtered_sam2"
rm "$in".amb "$in".ann "$in".bwt "$in".pac "$in".sa

printf "\n"
