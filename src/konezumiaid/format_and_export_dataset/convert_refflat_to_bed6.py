from __future__ import annotations
import pandas as pd
from pathlib import Path


def convert_refFlat_to_bed6(df_refflat: pd.DataFrame, output_file: Path) -> None:
    with open(output_file, "w") as outfile:
        for _, row in df_refflat.iterrows():
            chrom = row["chrom"]
            txStart = row["txStart"]
            txEnd = row["txEnd"]
            transcriot_name = row["name"]
            score = "0"
            strand = row["strand"]

            bed_line = "\t".join([chrom, txStart, txEnd, transcriot_name, score, strand])
            outfile.write(bed_line + "\n")
