from __future__ import annotations
import re


def search_c_to_t_guide_seq(orf_seq: str) -> list[str]:
    have_pam = re.findall(r"(?=(\w{21}GG))", orf_seq)
    # check "CAA", "CAG", or "CGA"'s "C" in 17-19bp positions from PAM
    editable_gRNAs = [
        match
        for match in have_pam
        if any(match[i : i + 3] in {"CAA", "CAG", "CGA"} for i in range(1, 4))
    ]
    return list(dict.fromkeys(editable_gRNAs))


def search_g_to_a_guide_seq(orf_seq: str) -> list[str]:
    have_pam = re.findall(r"(?=(CC\w{21}))", orf_seq)
    # check "TGG"'s "G" in 17-20bp positions from PAM
    editable_gRNAs = [
        match
        for match in have_pam
        if any(match[i : i + 3] == "TGG" for i in range(17, 21))
    ]
    return list(dict.fromkeys(editable_gRNAs))
