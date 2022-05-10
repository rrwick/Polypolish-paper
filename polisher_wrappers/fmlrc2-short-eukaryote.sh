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
echo "FMLRC2 POLISHING WITH SHORT READS (EUKARYOTE OPTIONS)"
echo "-----------------------------------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# Build the multi-string Burrows Wheeler transform.
bwt=fmlrc2_"$RANDOM".npy
gunzip -c "$r1" "$r2" | awk 'NR % 4 == 2' | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert "$bwt"

# Polish!
fmlrc2 -t 16 -k 21 59 80 -f 0 "$bwt" "$in" "$out".fasta

# Clean up.
rm "$bwt"

printf "\n"
