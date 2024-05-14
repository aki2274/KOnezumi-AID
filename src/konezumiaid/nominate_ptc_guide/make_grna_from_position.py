from __future__ import annotations
import re
from konezumiaid.create_gene_dataclass import GeneData


def extract_c_to_t_grna_from_position(
    transcript_record: GeneData, positions: list[int]
) -> list[str]:
    # Convert positions to grna sequence.
    ct_grna = []
    for editable_C_position in positions:
        # pro_sgRNA is 25bp sequence. Because the editable "C" is in 3 patterns. The posistions are -17~-19bp from PAM.
        pro_sgRNA = transcript_record.orf_seq[
            editable_C_position - 1 : editable_C_position + 24
        ]
        # extract the 20bp + PAM, as grna sequence
        rev_pro_sgRNA = pro_sgRNA[::-1]
        PAM_positions = list(re.finditer("(?=(GG))", rev_pro_sgRNA))
        for PAM_position in PAM_positions:
            correction_by_target_window = PAM_position.start()
            if 2 <= correction_by_target_window <= 4:
                sgRNA = transcript_record.orf_seq[
                    editable_C_position
                    + 1
                    - correction_by_target_window : editable_C_position
                    + 24
                    - correction_by_target_window
                ]
                # Check if the sgRNA has "TTTT" in the guide (first 20bp, excluding PAM).
                # If it has "TTTT", U6 promoter will not work.
                if "TTTT" not in sgRNA[:20]:
                    ct_grna.append(sgRNA)
    return ct_grna


def extract_g_to_a_grna_from_position(
    transcript_record: GeneData, positions: list[int]
) -> list[str]:
    #  convert positions to grna sequence.
    ga_grna = []
    for editable_codons_T_positon in positions:
        pro_sgRNA = transcript_record.orf_seq[
            editable_codons_T_positon - 20 : editable_codons_T_positon + 3
        ]
        PAM_positions = re.finditer("(?=(CC))", pro_sgRNA)
        # extract the grna sequence from the candidate sequence.
        for PAM_position in PAM_positions:
            correction_by_target_window = PAM_position.start()
            if 0 <= int(correction_by_target_window) <= 3:
                sgRNA = transcript_record.orf_seq[
                    editable_codons_T_positon
                    - 20
                    + correction_by_target_window : editable_codons_T_positon
                    + 3
                    + correction_by_target_window
                ]
                # Check if the sgRNA has "AAAA" in the guide (excluding PAM).
                # If it has "AAAA", U6 promoter will not work.
                if "AAAA" not in sgRNA[3:]:
                    ga_grna.append(sgRNA)
    return ga_grna
