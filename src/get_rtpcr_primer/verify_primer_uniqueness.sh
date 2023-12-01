#!/bin/bash

input_primers=/tmp/candidateprimer.
output_uniqueness=/mnt/c/KOnezumi-AID/data #ユーザーに合わせて変えるべき？

bowtie -a -v 0 -x data/indexes/mm39 -f data/reads/primers.fa > /mnt/c/KOnezumi-AID/data/boetie_result_0.sam
bowtie -a -v 1 -x data/indexes/mm39 -f data/reads/primers.fa > /mnt/c/KOnezumi-AID/data/boetie_result_1.sam
bowtie -a -v 2 -x data/indexes/mm39 -f data/reads/primers.fa > /mnt/c/KOnezumi-AID/data/boetie_result_2.sam

uniq_data_0="/mnt/c/KOnezumi-AID/data/boetie_result_0.sam"
uniq_data_1="/mnt/c/KOnezumi-AID/data/boetie_result_1.sam"
uniq_data_2="/mnt/c/KOnezumi-AID/data/boetie_result_2.sam"

awk '!/^@/{print $1}' "$uniq_data_0" | sort | uniq -c > output_sam_0.txt
awk '!/^@/{print $1}' "$uniq_data_1" | sort | uniq -c > output_sam_1.txt
awk '!/^@/{print $1}' "$uniq_data_2" | sort | uniq -c > output_sam_2.txt