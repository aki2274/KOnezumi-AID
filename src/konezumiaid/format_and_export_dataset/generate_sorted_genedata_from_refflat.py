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
    unique_df = df_refflat.drop_duplicates("name", keep=False)
    # remove the gene symbols, if the any of the transcript names are duplicated
    duplicates_df = df_refflat[df_refflat.duplicated("name", keep=False)]
    export_df = unique_df[~unique_df["geneName"].isin(duplicates_df["geneName"])]
    return export_df


def remove_NR_transcripts(df_refflat: pd.DataFrame) -> pd.DataFrame:
    # remove the transcripts that are not in the RefSeq database
    return df_refflat[~df_refflat["name"].str.contains("NR_")]


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


def convert_str_int_in_df(refflat_object: pd.DataFrame) -> pd.DataFrame:
    refflat_object["txStart"] = int(refflat_object["txStart"])
    refflat_object["txEnd"] = int(refflat_object["txEnd"])
    refflat_object["cdsStart"] = int(refflat_object["cdsStart"])
    refflat_object["cdsEnd"] = int(refflat_object["cdsEnd"])
    return refflat_object


def clea_by_txStart(refflta_object: pd.DataFrame) -> pd.DataFrame:
    """txStart is the start of the transcript.So set the start of the transcript to 0."""

    refflta_object["exonStarts"] = refflta_object["exonStarts"].split(",")
    refflta_object["exonEnds"] = refflta_object["exonEnds"].split(",")

    refflta_object = convert_str_int_in_df(refflta_object)

    refflta_object["txEnd"] = refflta_object["txEnd"] - refflta_object["txStart"]
    refflta_object["cdsStart"] = refflta_object["cdsStart"] - refflta_object["txStart"]
    refflta_object["cdsEnd"] = refflta_object["cdsEnd"] - refflta_object["txStart"]
    refflta_object["exonEnds"] = [
        int(element) - refflta_object["txStart"]
        for element in refflta_object["exonEnds"]
    ]
    refflta_object["exonStarts"] = [
        int(element) - refflta_object["txStart"]
        for element in refflta_object["exonStarts"]
    ]
    refflta_object["txStart"] = 0
    refflta_object["exonStarts"] = ",".join(
        map(str, sorted(refflta_object["exonStarts"]))
    )
    refflta_object["exonEnds"] = ",".join(map(str, sorted(refflta_object["exonEnds"])))
    return refflta_object


def clean_by_strand(refflat_object: pd.DataFrame) -> pd.DataFrame:
    """If the strand is "-", reverse the start and end of the transcript."""
    if refflat_object["strand"] == "-":
        refflat_object["exonStarts"] = refflat_object["exonStarts"].split(",")
        refflat_object["exonEnds"] = refflat_object["exonEnds"].split(",")
        refflat_object = convert_str_int_in_df(refflat_object)

        refflat_object["txStart"] = abs(
            refflat_object["txStart"] - refflat_object["txEnd"]
        )
        refflat_object["cdsStart"] = abs(
            refflat_object["cdsStart"] - refflat_object["txEnd"]
        )
        refflat_object["cdsEnd"] = abs(
            refflat_object["cdsEnd"] - refflat_object["txEnd"]
        )
        refflat_object["exonEnds"] = [
            abs(int(element) - refflat_object["txEnd"])
            for element in refflat_object["exonEnds"]
        ]
        refflat_object["exonStarts"] = [
            abs(int(element) - refflat_object["txEnd"])
            for element in refflat_object["exonStarts"]
        ]
        refflat_object["txEnd"] = 0

        txStart = refflat_object["txStart"]
        txEnd = refflat_object["txEnd"]
        cdsStart = refflat_object["cdsStart"]
        cdsEnd = refflat_object["cdsEnd"]
        exonStarts = refflat_object["exonStarts"]
        exonEnds = refflat_object["exonEnds"]

        refflat_object["txStart"] = txEnd
        refflat_object["txEnd"] = txStart
        refflat_object["cdsStart"] = cdsEnd
        refflat_object["cdsEnd"] = cdsStart
        refflat_object["exonStarts"] = exonEnds
        refflat_object["exonEnds"] = exonStarts
        refflat_object["exonStarts"] = ",".join(
            map(str, sorted(refflat_object["exonStarts"]))
        )
        refflat_object["exonEnds"] = ",".join(
            map(str, sorted(refflat_object["exonEnds"]))
        )
    return refflat_object


def clean_refflat(df_refflat: pd.DataFrame) -> pd.DataFrame:
    df_refflat_sorted = df_refflat.copy()
    df_refflat_sorted = df_refflat.apply(clea_by_txStart, axis=1)
    df_refflat_sorted = df_refflat_sorted.apply(clean_by_strand, axis=1)
    return df_refflat_sorted
