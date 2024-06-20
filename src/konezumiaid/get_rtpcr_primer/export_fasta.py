from __future__ import annotations
from pathlib import Path


def export_fasta(candidate_primers: list[dict], output_path: Path) -> None:
    """
    Export candidate primer sequences to a FASTA file for executing bowtie.
    Args:
        candidate_primers (list[dict]): A list of dictionaries containing candidate primer information.
        output_path (Path): The path to the output FASTA file.

    Returns:
        None

    Note; The output fasta file contains the same headers and sequences for each entry.

    Example:
        input:
        [{"left_seq": "AAAA", "right_seq": "TACG"},{"left_seq": "AAAA", "right_seq": "CATG"}]

        outputfile:
        >AAAA
        AAAA
        >TACG
        TACG
        >GCTA
        GCTA
    """
    processed_names = set()

    with open(output_path, "w") as file:
        for primer_pair in candidate_primers:
            left_primer = primer_pair["left"]
            right_primer = primer_pair["right"]

            left_primer_name = left_primer
            if left_primer_name not in processed_names:
                processed_names.add(left_primer_name)
                file.write(f">{left_primer_name}\n{left_primer}\n")
            right_primer_name = right_primer
            if right_primer_name not in processed_names:
                processed_names.add(right_primer_name)
                file.write(f">{right_primer_name}\n{right_primer}\n")
