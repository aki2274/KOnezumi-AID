import re
import pandas as pd


def read_fasta(path: str) -> str:
    sequences = {}
    seq_id = None
    current_seq = ""
    with open(path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if seq_id is not None:
                    sequences[seq_id] = current_seq
                    current_seq = ""
                seq_id = line[1:]
            else:
                current_seq += line
        if seq_id is not None and current_seq != "":
            sequences[seq_id] = current_seq
    return sequences


def read_refFlat(path: str) -> str:
    with open(path, "r", encoding="utf-8") as file:
        lines = file.readlines()
    result = []
    for line in lines:
        line = re.split("\s+", line)
        contents = [line[i] for i in range(len(line) - 1)]
        result.append(contents)
    return result


def make_dataset(path_refFlat):
    gene_data = read_refFlat(path_refFlat)
    col_names = [
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
    gene_df = pd.DataFrame(gene_data, columns=col_names)
    return gene_df
