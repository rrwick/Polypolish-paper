Polishing tools often require multiple commands to run (e.g. aligner indexing, read alignment and then polishing), so the Bash scripts in this directory wrap polishing commands in in an easy and consistent manner. These scripts take of the following (as necessary):
* make sure the input files exists
* create an aligner index file
* align the reads
* index the alignment file
* run the polisher
* convert the output to uppercase and one-line-per-sequence FASTA format
* delete intermediate files

Short-read polishing scripts are named `*-short.sh` and are run like this:
```bash
./tool-short.sh input.fasta reads_1.fastq.gz reads_2.fastq.gz output_prefix
```

Hybrid polishing scripts are named `*-hybrid.sh` and are run like this:
```bash
./tool-hybrid.sh input.fasta reads_1.fastq.gz reads_2.fastq.gz reads_long.fastq.gz output_prefix
```

The result is a file named `output_prefix.fasta` which contains the polished sequence.
