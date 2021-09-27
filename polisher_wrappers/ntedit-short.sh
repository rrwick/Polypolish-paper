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
echo "NTEDIT POLISHING WITH SHORT READS"
echo "---------------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# Polish!
echo -e "$r1\n$r2" > reads.in
nthits -b 36 -k 40 -t 16 --outbloom --solid @reads.in
ntedit -f "$in" -r solids_k40.bf -t 16 -m 1
mv *_rsolids_k40.bf_*_edited.fa "$out".fasta

# Make uppercase and one line per sequence.
temp=temp_"$RANDOM".fasta
seqtk seq "$out".fasta | awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' > "$temp"
mv "$temp" "$out".fasta

# Clean up.
rm *_rsolids_*.vcf *_rsolids_*.tsv
rm reads.in
rm solids_k40.bf

printf "\n"
