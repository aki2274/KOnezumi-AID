from __future__ import annotations


def create_list(cand_grna: list[str]) -> list[dict]:
    grna_list = []
    for grna in cand_grna:
        grna_list.append({"seq": grna})
    return grna_list
