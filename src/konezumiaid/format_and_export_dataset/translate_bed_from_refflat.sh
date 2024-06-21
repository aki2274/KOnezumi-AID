#!/bin/bash
# Execute bedtools getfasta command
bedtools getfasta -nameOnly -fi $1 -bed $2 -fo $3