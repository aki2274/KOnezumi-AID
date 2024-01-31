#!/bin/bash

# Execute bedtools getfasta command
bedtools getfasta -fi $1 -bed $2 -name+
