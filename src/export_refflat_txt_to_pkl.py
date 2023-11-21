from __future__ import annotations
import pandas as pd
import pickle


def read_refFlat(path_refFlat: str) -> list[str]:
    result = []
    with open(path_refFlat, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            # Remove the trailing comma
            line = [element.strip(",") for element in line]
            result.append(line)
    return result


def make_genedata_dataframe(path_refFlat: str) -> pd.DataFrame:
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
    return df


def export_genedata_pkl(
    gene_dataframe: pd.DataFrame, path_output: str = "refFlat_genedata.pkl"
) -> None:
    # Remove duplicates of transcriptions
    duplicates_df = gene_dataframe[gene_dataframe.duplicated("name", keep=False)]
    unique_df = gene_dataframe.drop_duplicates("name", keep=False)
    export_df = unique_df[~unique_df["geneName"].isin(duplicates_df["geneName"])]
    export_data = export_df.to_dict(orient="records")
    with open(path_output, "wb") as file:
        pickle.dump(export_data, file)