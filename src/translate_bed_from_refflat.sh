#!/bin/bash
# Execute bedtools getfasta command
bedtools getfasta -name+ -fi $1 -bed $2 -fo $3
