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



# This is the BWA-single-end-read-based version of Racon short-read polishing that I tried but ultimately didn't use:

# # Prepare a single file of short reads. BWA drops the  '/1' and '/2' at the end of read names,
# # which upsets Racon (because the SAM read names don't match the FASTQ read names), so I replace
# # those with '_1' and '_2'
# gunzip -c "$r1" | sed 's|/1$|_1|' > combined_reads.fastq
# gunzip -c "$r2" | sed 's|/1$|_2|' >> combined_reads.fastq

# # Prepare the short-read alignments.
# short_alignments="$in".sam
# bwa index "$in"
# bwa mem -t 16 "$in" combined_reads.fastq > "$short_alignments"

# # Polish!
# racon -t 16 --no-trimming combined_reads.fastq "$short_alignments" "$in" > "$out".fasta

# # Clean up.
# rm combined_reads.fastq
# rm "$short_alignments"
# rm "$in".amb "$in".ann "$in".bwt "$in".pac "$in".sa




# This is the BWA-paired-end-read-based version of Racon short-read polishing that I tried but ultimately didn't use:

# # Prepare a single file of short reads.
# gunzip -c "$r1" > combined_reads.fastq
# gunzip -c "$r2" >> combined_reads.fastq

# # Prepare the short-read alignments, using BWA in a paired-read manner.
# short_alignments="$in".sam
# bwa index "$in"
# bwa mem -t 16 "$in" "$r1" "$r2" > "$short_alignments".temp

# # BWA removes '/1' and '/2' from read names, and Racon doesn't like this: it thinks the read names
# # in the SAM aren't in the FASTQ. So I use some samtools/awk to add '/1' to the first-in-pair reads
# # and '/2' to the second-in-pair reads.
# samtools view -H "$short_alignments".temp > "$short_alignments"
# samtools view -h -f 64 "$short_alignments".temp | awk -F $'\t' 'BEGIN {OFS = FS} {$1=$1"/1"}1' >> "$short_alignments"
# samtools view -f 128 "$short_alignments".temp | awk -F $'\t' 'BEGIN {OFS = FS} {$1=$1"/2"}1' >> "$short_alignments"

# # Polish!
# racon -t 16 --no-trimming combined_reads.fastq "$short_alignments" "$in" > "$out".fasta

# # Clean up.
# rm combined_reads.fastq
# rm "$short_alignments".temp "$short_alignments"
# rm "$in".amb "$in".ann "$in".bwt "$in".pac "$in".sa
