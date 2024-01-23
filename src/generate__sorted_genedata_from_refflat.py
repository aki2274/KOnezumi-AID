from __future__ import annotations
import pandas as pd
import pickle


def read_refFlat(path_refFlat: str) -> list[list[str]]:
    # Read the refFlat file and return a list of lists.
    result = []
    with open(path_refFlat, "r") as file:
        for line in file:
            # Split the line and strip trailing commas from each element
            line = [element.strip(",") for element in line.strip().split("\t")]
            result.append(line)
    return result


def built_gene_dataframe(path_refFlat: str) -> pd.DataFrame:
    # label the columns according to the UCSC refFlat file
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
    return pd.DataFrame(results)


def sort_transcript_data(genedata: pd.DataFrame) -> pd.DataFrame:
    # sort the transcript data according to the txStart
    sorted_data = genedata.copy()
    sorted_data["exonStarts"] = sorted_data["exonStarts"].split(",")
    sorted_data["exonEnds"] = sorted_data["exonEnds"].split(",")

    sorted_data["txEnd"] = int(sorted_data["txEnd"]) - int(sorted_data["txStart"])
    sorted_data["cdsStart"] = int(sorted_data["cdsStart"]) - int(sorted_data["txStart"])
    sorted_data["cdsEnd"] = int(sorted_data["cdsEnd"]) - int(sorted_data["txStart"])
    sorted_data["exonEnds"] = [
        int(element) - int(sorted_data["txStart"])
        for element in sorted_data["exonEnds"]
    ]
    sorted_data["exonStarts"] = [
        int(element) - int(sorted_data["txStart"])
        for element in sorted_data["exonStarts"]
    ]
    sorted_data["txStart"] = int(sorted_data["txStart"]) - int(sorted_data["txStart"])
    if sorted_data["strand"] == "-":
        # if the transcript is on the negative strand, reverse start and end.
        sorted_data["txStart"] = abs(sorted_data["txStart"] - sorted_data["txEnd"])
        sorted_data["cdsStart"] = abs(sorted_data["cdsStart"] - sorted_data["txEnd"])
        sorted_data["cdsEnd"] = abs(sorted_data["cdsEnd"] - sorted_data["txEnd"])
        sorted_data["exonEnds"] = [
            abs(element - sorted_data["txEnd"]) for element in sorted_data["exonEnds"]
        ]
        sorted_data["exonStarts"] = [
            abs(element - sorted_data["txEnd"]) for element in sorted_data["exonStarts"]
        ]
        sorted_data["txEnd"] = sorted_data["txEnd"] - sorted_data["txEnd"]

        txStart = sorted_data["txStart"]
        txEnd = sorted_data["txEnd"]
        cdsStart = sorted_data["cdsStart"]
        cdsEnd = sorted_data["cdsEnd"]
        exonStarts = sorted_data["exonStarts"]
        exonEnds = sorted_data["exonEnds"]

        sorted_data["txStart"] = txEnd
        sorted_data["txEnd"] = txStart
        sorted_data["cdsStart"] = cdsEnd
        sorted_data["cdsEnd"] = cdsStart
        sorted_data["exonStarts"] = exonEnds
        sorted_data["exonEnds"] = exonStarts
    sorted_data["exonStarts"] = ",".join(map(str, sorted(sorted_data["exonStarts"])))
    sorted_data["exonEnds"] = ",".join(map(str, sorted(sorted_data["exonEnds"])))
    return sorted_data


def sort_gene_dataframe(gene_dataframe: pd.DataFrame) -> pd.DataFrame:
    gene_dataframe = gene_dataframe.apply(sort_transcript_data, axis=1)
    return gene_dataframe


def remove_genename_duplicates(
    gene_dataframe: pd.DataFrame,
) -> None:
    # remove the duplicates of transcript names
    duplicates_df = gene_dataframe[gene_dataframe.duplicated("name", keep=False)]
    unique_df = gene_dataframe.drop_duplicates("name", keep=False)
    export_df = unique_df[~unique_df["geneName"].isin(duplicates_df["geneName"])]
    return export_df
