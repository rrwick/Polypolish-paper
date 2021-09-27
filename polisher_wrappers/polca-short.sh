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
echo "POLCA POLISHING WITH SHORT READS"
echo "--------------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# Polish!
/home/ubuntu/programs/MaSuRCA-4.0.3/bin/polca.sh -a "$in" -r $r1" "$r2 -t 16 -m 1G
mv *.PolcaCorrected.fa "$out".fasta

# Make uppercase and one line per sequence.
temp=temp_"$RANDOM".fasta
seqtk seq "$out".fasta | awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' > "$temp"
mv "$temp" "$out".fasta

# Clean up.
rm "$in".fai
rm *.err *.bwa.* *.sam *.success *.bam *.bam.bai *.batches *.names *.report *.vcf

printf "\n"
