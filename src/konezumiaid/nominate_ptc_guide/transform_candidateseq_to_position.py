from __future__ import annotations
import re


def transform_c_to_t_guideseq_to_position(
    orf_seq: str, guide_seq: list[str]
) -> list[int]:
    # Get the position of "CAA", "CAG", or "CGA"'s "C" in target-window.(the position on the ORF)
    editable_c_positions = []
    for seq in guide_seq:
        for seq_position in re.finditer(seq, orf_seq):
            # target gRNA is 20bp, but target window 17-19pb upstream from PAM. So use num "1"
            position = (
                re.search(r"(CAA|CAG|CGA)", seq[1:]).start() + seq_position.start() + 1
            )
            editable_c_positions.append(position)
    editable_c_positions = list(dict.fromkeys(editable_c_positions))
    return editable_c_positions


def transform_g_to_a_guideseq_to_position(
    orf_seq: str, guide_seq: list[str]
) -> list[int]:
    # Get the position of "TGG"'s "T", if "TGG"'s "G" in target-window.(the position on the ORF)
    editable_codon_T_positions = []
    for seq in guide_seq:
        for seq_position in re.finditer(seq, orf_seq):
            rev_seq = seq[::-1]
            if rev_seq[0:3] == "GGT" and rev_seq[3:6] == "GGT":
                # The case that guide-sequence CCN NNNNNNNNNNNNNNTGGTGG
                editable_codon_T_positions.append(seq_position.start() + len(seq) - 3)
                editable_codon_T_positions.append(
                    seq_position.start() + len(seq) - 3 - 3
                )
            else:
                correction_by_target_window = re.search("GGT", rev_seq).start()
                editable_codon_T_positions.append(
                    seq_position.start() + len(seq) - 3 - correction_by_target_window
                )
    editable_codon_T_positions = list(dict.fromkeys(editable_codon_T_positions))
    return editable_codon_T_positions
