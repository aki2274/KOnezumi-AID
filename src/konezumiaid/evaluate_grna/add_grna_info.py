from __future__ import annotations


def link_to_crisperdirect(candidates: list[dict]) -> list[dict]:
    for candidate in candidates:
        candidate["link_to_crisprdirect"] = f"https://crispr.dbcls.jp/?userseq={candidate['seq']}&pam=NGG&db=mm39"
    return candidates
