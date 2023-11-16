from __future__ import annotations
import pandas as pd


def read_refFlat(path_refFlat: str) -> list[str]:
    result = []
    with open(path_refFlat, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            # Remove the trailing comma
            line = [element.strip(",") for element in line]
            result.append(line)
    return result


def export_genedata_pkl(
    path_refFlat: str, path_output: str = "refFlat_genedata.pkl"
) -> None:
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
    export_data.to_pickle(path_output)
