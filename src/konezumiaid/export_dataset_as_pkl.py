from __future__ import annotations
from pathlib import Path
import pickle
import subprocess
from konezumiaid.convert_refflat_to_bed import convert_refFlat_to_bed
from konezumiaid.generate_seq_dict_from_fasta import (
    read_fasta,
    create_sorted_seq_dict,
)
from konezumiaid.generate_sorted_genedata_from_refflat import (
    built_gene_dataframe,
    sort_gene_dataframe,
    remove_genename_duplicates,
)


Path("data").mkdir(parents=True, exist_ok=True)


def export_pkl(refflat_path: Path, fasta_path: Path) -> None:
    """
    Export the refflat file and the sequence dictionary as pickle files.

    Args:
        refflat_path (Path): path to the refflat.txt file.
        fasta_path (Path): path to the mm39 fasta file.
        out_refflat_path (Path): path to the output refflat file.
        out_dict_path (Path): path to the output sequence dictionary file.

    Returns:
        None
    example:
        >>> refflat_path = Path("data", "refFlat.txt")
        >>> fasta_path = Path("data", "mm39.fa")
        >>> out_refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
        >>> out_dict_path = Path("data", "sorted_seq_dict.pkl")
        >>> export_pkl(refflat_path, fasta_path, out_refflat_path, out_dict_path)
    then, the sorted files are exported as pickle.
    """
    out_refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    out_dict_path = Path("data", "sorted_seq_dict.pkl")

    bed_output_path = Path("data", "refFlat.bed")
    bed_fast_path = Path("data", "bed_refFlat.fa")
    convert_refFlat_to_bed(refflat_path, bed_output_path)
    translate_bed_path = Path("src", "translate_bed_from_refflat.sh")
    subprocess.run(
        [
            "bash",
            translate_bed_path,
            str(fasta_path),
            str(bed_output_path),
            str(bed_fast_path),
        ]
    )
    gene_data = built_gene_dataframe(refflat_path)
    gene_data = remove_genename_duplicates(gene_data)
    sorted_gene_data = remove_genename_duplicates(sort_gene_dataframe(gene_data))
    gene_seq = read_fasta(bed_fast_path)
    seq_dict = create_sorted_seq_dict(gene_data, sorted_gene_data, gene_seq)
    with open(out_dict_path, "wb") as f:
        pickle.dump(seq_dict, f)
    refflat_data = sorted_gene_data.to_dict(orient="records")
    with open(out_refflat_path, "wb") as f:
        pickle.dump(refflat_data, f)
