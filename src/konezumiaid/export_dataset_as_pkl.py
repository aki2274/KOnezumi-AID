from __future__ import annotations
from pathlib import Path
import pickle
import subprocess
from konezumiaid.convert_refflat_to_bed6 import convert_refFlat_to_bed6
from konezumiaid.generate_seq_dict_from_fasta import (
    read_fasta,
    create_strand_plus_seq_dict,
)
from konezumiaid.generate_sorted_genedata_from_refflat import (
    built_gene_dataframe,
    sort_gene_dataframe,
    remove_transcript_duplicates,
)


Path("data").mkdir(parents=True, exist_ok=True)


def export_pkl(refflat_path: Path, chromosome_fasta_path: Path) -> None:
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
    sorted_refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    sorted_transcript_seq_dict_path = Path("data", "sorted_seq_dict.pkl")

    bed6_path = Path("data", "refFlat.bed")
    convert_refFlat_to_bed6(refflat_path, bed6_path)

    transcripts_fast_path = Path("data", "bed_refFlat.fa")
    bedtools_path = Path("src", "konezumiaid", "translate_bed_from_refflat.sh")
    subprocess.run(
        [
            "bash",
            bedtools_path,
            str(chromosome_fasta_path),
            str(bed6_path),
            str(transcripts_fast_path),
        ]
    )
    df_refflat = built_gene_dataframe(refflat_path)
    df_refflat = remove_transcript_duplicates(df_refflat)
    df_refflat_sorted = remove_transcript_duplicates(sort_gene_dataframe(df_refflat))
    transcript_seq_dict = read_fasta(transcripts_fast_path)
    transcript_seq_sorted_dict = create_strand_plus_seq_dict(
        df_refflat, df_refflat_sorted, transcript_seq_dict
    )
    with open(sorted_transcript_seq_dict_path, "wb") as f:
        pickle.dump(transcript_seq_sorted_dict, f)
    sorted_refflat = df_refflat_sorted.to_dict(orient="records")
    with open(sorted_refflat_path, "wb") as f:
        pickle.dump(sorted_refflat, f)
