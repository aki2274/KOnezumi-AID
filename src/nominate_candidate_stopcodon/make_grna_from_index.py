from __future__ import annotations
from dataclasses import dataclass
import re


@dataclass
class GeneData:
    orf_seq: str
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


def convert_ct_grna(ds: GeneData, indices: list[int]) -> list[str]:
    # Convert indices to grna sequence.
    targets = [index - ds.txStart for index in indices]
    ct_grna = []
    # I don't know where the target "C" is in -17~-19bp. So, I get the longer canditdae sequence.
    for target in targets:
        pro_sgRNA = ds.orf_seq[target - 1 : target + 24]
        reversed_pro_sgRNA = pro_sgRNA[::-1]
        matches = list(re.finditer("(?=(GG))", reversed_pro_sgRNA))
        # extract the grna sequence from the candidate sequence.
        for match in matches:
            start_position = match.start()
            if 2 <= start_position <= 4:
                sgRNA = ds.orf_seq[
                    target + 1 - start_position : target + 24 - start_position
                ]
                ct_grna.append(sgRNA)
    return ct_grna


def convert_ga_grna(ds: GeneData, indices: list[int]) -> list[str]:
    #  convert index to grna sequence.
    targets = [index - ds.txStart for index in indices]
    ga_grna = []
    # I don't know where the target "G" is in 17~19bp. So, I get the longer candidate sequence.
    for target in targets:
        pro_sgRNA = ds.orf_seq[target - 20 : target + 3]
        matches = re.finditer("(?=(CC))", pro_sgRNA)
        # extract the grna sequence from the candidate sequence.
        for match in matches:
            add_num = match.start()
            if 0 <= int(add_num) <= 3:
                sgRNA = ds.orf_seq[target - 20 + add_num : target + 3 + add_num]
                ga_grna.append(sgRNA)
    return ga_grna
