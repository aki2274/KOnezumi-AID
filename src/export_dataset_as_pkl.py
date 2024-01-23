from __future__ import annotations
from pandas import DataFrame
from pathlib import Path
from src.generate_seq_dict_from_fasta import (
    read_fasta,
    create_sorted_seq_dict,
)
from src.generate__sorted_genedata_from_refflat import (
    built_gene_dataframe,
    sort_gene_dataframe,
    remove_genename_duplicates,
)
import pickle


def export_pkl(
    refflat_path: Path, fasta_path: Path, out_refflat_path: Path, out_dict_path
) -> None:
    gene_data = built_gene_dataframe(refflat_path)
    gene_data = remove_genename_duplicates(gene_data)
    sorted_gene_data = remove_genename_duplicates(sort_gene_dataframe(gene_data))
    gene_seq = read_fasta(fasta_path)
    seq_dict = create_sorted_seq_dict(gene_data, sorted_gene_data, gene_seq)
    with open(out_dict_path, "wb") as f:
        pickle.dump(seq_dict, f)
    refflat_data = gene_data.to_dict(orient="records")
    with open(out_refflat_path, "wb") as f:
        pickle.dump(refflat_data, f)
