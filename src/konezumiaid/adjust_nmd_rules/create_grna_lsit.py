from __future__ import annotations


def create_list(ct_cand: list[str], ga_cand: list[str]) -> list[dict]:
    grna_list = []
    for grna in ct_cand:
        grna_list.append({"ct_seq": grna})
    for grna in ga_cand:
        grna_list.append({"ga_seq": grna})
    return grna_list
