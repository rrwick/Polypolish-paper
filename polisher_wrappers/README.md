Polishing tools often require multiple commands to run (e.g. aligner indexing, read alignment and then polishing), so the Bash scripts in this directory wrap polishing commands in in an easy and consistent manner.

All short-read polishing scripts are named `*-short.sh` and are run like this:
```bash
tool-short.sh input.fasta reads_1.fastq.gz reads_2.fastq.gz output_prefix
```

All hybrid polishing scripts are named `*-hybrid.sh` and are run like this:
```bash
tool-short.sh input.fasta reads_1.fastq.gz reads_2.fastq.gz reads_long.fastq.gz output_prefix
```

These scripts take of the following (as necessary):
* make sure the input files exists
* create an aligner index file
* align the reads
* index the alignment file
* run the polisher
* convert the output to uppercase and one-line-per-sequence FASTA format
* delete intermediate files
