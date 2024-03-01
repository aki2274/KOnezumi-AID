#!/bin/bash
input_primers="data/uniq/candidateprimer.fa"
bowtie -a -v 0 -x data/indexes/mm39 -f "$input_primers"  > data/uniq/bowtie_result_0.sam
bowtie -a -v 1 -x data/indexes/mm39 -f "$input_primers"  > data/uniq/bowtie_result_1.sam
bowtie -a -v 2 -x data/indexes/mm39 -f "$input_primers"  > data/uniq/bowtie_result_2.sam
grep -v "^@" data/uniq/bowtie_result_0.sam | cut -f 1 | sort | uniq -c > data/uniq/0_miss_counts.txt
grep -v "^@" data/uniq/bowtie_result_1.sam | cut -f 1 | sort | uniq -c > data/uniq/1_miss_counts.txt
grep -v "^@" data/uniq/bowtie_result_2.sam | cut -f 1 | sort | uniq -c > data/uniq/2_miss_counts.txt
