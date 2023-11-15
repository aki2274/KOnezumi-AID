from __future__ import annotations
import re
import pandas as pd


def read_refFlat(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as file:
        lines = file.readlines()
    result = []
    for line in lines:
        line = re.split("\s+", line)
        contents = [element.rstrip(",") for element in line[:-1]]
        result.append(contents)
    return result


def export_genedata_csv(path_refFlat: str):
    values_refflat = read_refFlat(path_refFlat)
    key_names = [
        "geneName",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
    ]
    results = [dict(zip(key_names, value)) for value in values_refflat]
    df = pd.DataFrame(results)
    # Remove duplicates of transcriptions
    duplicates_df = df[df.duplicated("name", keep=False)]
    unique_df = df.drop_duplicates("name", keep=False)
    export_data = unique_df[~unique_df["geneName"].isin(duplicates_df["geneName"])]
    export_data.to_csv("refFlat_genedata.csv", index=False)
