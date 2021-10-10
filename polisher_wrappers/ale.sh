#!/usr/bin/env bash

in=$1   # FASTA to assess
r1=$2   # short reads (first in pairs)
r2=$3   # short reads (second in pairs)
out=$4  # output prefix

# Make sure the necessary files exist.
if [[ ! -f "$in" ]] ; then echo $in" does not exist"; exit 1; fi
if [[ ! -f "$r1" ]] ; then echo $r1" does not exist"; exit 1; fi
if [[ ! -f "$r2" ]] ; then echo $r2" does not exist"; exit 1; fi

printf "\n"
echo "ASSESSING ASSEMBLY WITH ALE"
echo "---------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# Prepare the short-read alignments.
sam=ale_"$RANDOM".sam
bwa index "$in"
bwa mem -t 16 "$in" "$r1" "$r2" > "$sam"

# Assess!
ALE "$sam" "$in" "$out".ale

# Clean up.
rm "$sam"
rm "$in".amb "$in".ann "$in".bwt "$in".pac "$in".sa
rm "$out".ale.param

printf "\n"

