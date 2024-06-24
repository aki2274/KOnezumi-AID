from __future__ import annotations
from pathlib import Path
import pickle
import subprocess
from konezumiaid.format_and_export_dataset.convert_refflat_to_bed6 import (
    convert_refFlat_to_bed6,
)
from konezumiaid.format_and_export_dataset.generate_seq_dict_from_fasta import (
    read_fasta,
    create_strand_plus_seq_dict,
)
from konezumiaid.format_and_export_dataset.generate_sorted_genedata_from_refflat import (
    built_gene_dataframe,
    clean_refflat,
    remove_transcript_duplicates,
    remove_NR_transcripts,
)


def export_pkl(refflat_path: Path, chromosome_fasta_path: Path) -> None:
    """
    Format and export the dataset as pickle files.
    The starting point of the transcript is set to 0. remove NR transcripts and duplicates.
    Export the sorted refflat and sequence dictionary as pickle files.

    Args:
        refflat_path (Path): Path to the refFlat txt file.
        chromosome_fasta_path (Path): Path to the chromosome fasta file.(ex. mm39.fa)

    Returns:
        None
    """
    sorted_refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    sorted_transcript_seq_dict_path = Path("data", "sorted_seq_dict.pkl")

    bed6_path = Path("data", "refFlat.bed")

    df_refflat = built_gene_dataframe(refflat_path)
    df_refflat = remove_transcript_duplicates(df_refflat)
    df_removed = remove_NR_transcripts(df_refflat)
    convert_refFlat_to_bed6(df_removed, bed6_path)

    transcripts_fast_path = Path("data", "bed_refFlat.fa")
    bedtools_path = Path(
        "src",
        "konezumiaid",
        "format_and_export_dataset",
        "translate_bed_from_refflat.sh",
    )
    subprocess.run(
        [
            "bash",
            bedtools_path,
            str(chromosome_fasta_path),
            str(bed6_path),
            str(transcripts_fast_path),
        ]
    )
    print("converted bed to fasta.")
    df_refflat_sorted = df_removed.copy()
    df_refflat_sorted = clean_refflat(df_refflat_sorted)

    transcript_seq_dict = read_fasta(transcripts_fast_path)
    print("creating transcript_seq_dict")
    transcript_seq_sorted_dict = create_strand_plus_seq_dict(df_removed, transcript_seq_dict)
    print("created transcript_seq_dict")
    with open(sorted_transcript_seq_dict_path, "wb") as f:
        pickle.dump(transcript_seq_sorted_dict, f)
    sorted_refflat = df_refflat_sorted.to_dict(orient="records")
    with open(sorted_refflat_path, "wb") as f:
        pickle.dump(sorted_refflat, f)


def execute_export(refflat_path: Path, fasta_path: Path):
    if not refflat_path.exists() or not fasta_path.exists():
        raise FileNotFoundError("One or both of the specified files were not found.")
    if refflat_path.suffix != ".txt":
        raise ValueError("The refflat file must be .txt file.")
    if fasta_path.suffix != ".fa" or fasta_path.suffix != ".fasta":
        raise ValueError("The fasta file must be .fa or .fasta file.")
    export_pkl(refflat_path, fasta_path)
    print("Exported the dataset as pickle files.")
