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
echo "NEXTPOLISH POLISHING WITH SHORT READS"
echo "-------------------------------------"
echo "Input FASTA: "$in
echo "Input reads: "$r1", "$r2
printf "\n"

# First round of polishing (task 1)
bam=nextpolish_"$RANDOM".bam
bwa index "$in"
bwa mem -t 16 "$in" "$r1" "$r2" | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - "$bam"
samtools index -@ 16 "$bam"
samtools faidx "$in"
python /home/ubuntu/programs/NextPolish/lib/nextpolish1.py -g "$in" -t 1 -p 16 -s "$bam" -ploidy 1 > genome.polishtemp.fa

# Clean up
rm "$in".amb "$in".ann "$in".bwt "$in".pac "$in".sa "$in".fai
rm "$bam" "$bam".bai

# Second round of polishing (task 2)
in2=genome.polishtemp.fa
bwa index "$in2"
bwa mem -t 16 "$in2" "$r1" "$r2" | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - "$bam"
samtools index -@ 16 "$bam"
samtools faidx "$in2"
python /home/ubuntu/programs/NextPolish/lib/nextpolish1.py -g "$in2" -t 2 -p 16 -s "$bam" -ploidy 1 > "$out".fasta

# Make uppercase and one line per sequence.
temp=temp_"$RANDOM".fasta
seqtk seq "$out".fasta | awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' > "$temp"
mv "$temp" "$out".fasta

# Clean up
rm genome.polishtemp.fa
rm "$in2".amb "$in2".ann "$in2".bwt "$in2".pac "$in2".sa "$in2".fai
rm "$bam" "$bam".bai

printf "\n"
