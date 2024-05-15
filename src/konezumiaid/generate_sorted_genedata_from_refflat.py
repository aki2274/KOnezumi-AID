from __future__ import annotations
import pandas as pd


def read_refflat(path_refFlat: str) -> list[list[str]]:
    result = []
    with open(path_refFlat, "r") as file:
        for line in file:
            line = [element.strip(",") for element in line.strip().split("\t")]
            result.append(line)
    return result


def remove_transcript_duplicates(
    df_refflat: pd.DataFrame,
) -> pd.DataFrame:
    # remove the duplicates of transcript names
    duplicates_df = df_refflat[df_refflat.duplicated("name", keep=False)]
    unique_df = df_refflat.drop_duplicates("name", keep=False)
    export_df = unique_df[~unique_df["geneName"].isin(duplicates_df["geneName"])]
    return export_df


def built_gene_dataframe(refflat_path: str) -> pd.DataFrame:
    # label the columns according to the UCSC refFlat file
    values_refflat = read_refflat(refflat_path)
    column_names = [
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
    refflat_dict = [dict(zip(column_names, value)) for value in values_refflat]
    return pd.DataFrame(refflat_dict)


def convert_str_int_in_df(df_refflat: pd.DataFrame) -> pd.DataFrame:
    df_refflat["txStart"] = df_refflat["txStart"].astype(int)
    df_refflat["txEnd"] = df_refflat["txEnd"].astype(int)
    df_refflat["cdsStart"] = df_refflat["cdsStart"].astype(int)
    df_refflat["cdsEnd"] = df_refflat["cdsEnd"].astype(int)
    df_refflat["exonCount"] = df_refflat["exonCount"].astype(int)
    return df_refflat


def clea_by_txStart(df_refflat: pd.DataFrame) -> pd.DataFrame:
    """txStart is the start of the transcript.So set the start of the transcript to 0."""

    df_refflat["exonStarts"] = df_refflat["exonStarts"].split(",")
    df_refflat["exonEnds"] = df_refflat["exonEnds"].split(",")

    df_refflat = convert_str_int_in_df(df_refflat)

    df_refflat["txEnd"] = df_refflat["txEnd"] - df_refflat["txStart"]
    df_refflat["cdsStart"] = df_refflat["cdsStart"] - df_refflat["txStart"]
    df_refflat["cdsEnd"] = df_refflat["cdsEnd"] - df_refflat["txStart"]
    df_refflat["exonEnds"] = [
        int(element) - df_refflat["txStart"] for element in df_refflat["exonEnds"]
    ]
    df_refflat["exonStarts"] = [
        int(element) - df_refflat["txStart"] for element in df_refflat["exonStarts"]
    ]
    df_refflat["txStart"] = 0
    df_refflat["exonStarts"] = ",".join(map(str, sorted(df_refflat["exonStarts"])))
    df_refflat["exonEnds"] = ",".join(map(str, sorted(df_refflat["exonEnds"])))
    return df_refflat


def clean_by_strand(df_refflat: pd.DataFrame) -> pd.DataFrame:
    """If the strand is "-", reverse the start and end of the transcript."""

    df_refflat["exonStarts"] = df_refflat["exonStarts"].split(",")
    df_refflat["exonEnds"] = df_refflat["exonEnds"].split(",")

    if df_refflat["strand"] == "-":
        df_refflat["txStart"] = abs(df_refflat["txStart"] - df_refflat["txEnd"])
        df_refflat["cdsStart"] = abs(df_refflat["cdsStart"] - df_refflat["txEnd"])
        df_refflat["cdsEnd"] = abs(df_refflat["cdsEnd"] - df_refflat["txEnd"])
        df_refflat["exonEnds"] = [
            abs(int(element) - df_refflat["txEnd"])
            for element in df_refflat["exonEnds"]
        ]
        df_refflat["exonStarts"] = [
            abs(int(element) - df_refflat["txEnd"])
            for element in df_refflat["exonStarts"]
        ]
        df_refflat["txEnd"] = 0

        txStart = df_refflat["txStart"]
        txEnd = df_refflat["txEnd"]
        cdsStart = df_refflat["cdsStart"]
        cdsEnd = df_refflat["cdsEnd"]
        exonStarts = df_refflat["exonStarts"]
        exonEnds = df_refflat["exonEnds"]

        df_refflat["txStart"] = txEnd
        df_refflat["txEnd"] = txStart
        df_refflat["cdsStart"] = cdsEnd
        df_refflat["cdsEnd"] = cdsStart
        df_refflat["exonStarts"] = exonEnds
        df_refflat["exonEnds"] = exonStarts
    df_refflat["exonStarts"] = ",".join(map(str, sorted(df_refflat["exonStarts"])))
    df_refflat["exonEnds"] = ",".join(map(str, sorted(df_refflat["exonEnds"])))
    return df_refflat


def sort_gene_dataframe(df_refflat: pd.DataFrame) -> pd.DataFrame:
    df_refflat_sorted = df_refflat.apply(clea_by_txStart, axis=1)
    df_refflat_sorted = df_refflat_sorted.apply(clean_by_strand, axis=1)
    return df_refflat_sorted
