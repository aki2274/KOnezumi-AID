from __future__ import annotations
from pathlib import Path


def convert_refFlat_to_bed(refFlat_file: Path, output_file: Path) -> None:
    with open(refFlat_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            fields = line.strip().split("\t")

            chrom = fields[2]
            txStart = fields[4]
            txEnd = fields[5]
            name = fields[1]
            score = "0"
            strand = fields[3]

            bed_line = "\t".join([chrom, txStart, txEnd, name, score, strand])
            outfile.write(bed_line + "\n")
