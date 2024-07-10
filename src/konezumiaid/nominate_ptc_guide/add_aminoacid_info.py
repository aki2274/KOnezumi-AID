from __future__ import annotations


def transrate_codon_to_aminoacid(seq: str) -> str:
    """
    Translate the codon to amino acid.
    """
    codon_dict = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }
    seq = seq.upper()
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    amino_acid = ""
    for codon in codons:
        amino_acid += codon_dict[codon]
    return amino_acid


def link_position_and_aminoacid(positions_orf: list[int], position_cds: list[int], seq: list[str]) -> list[dict]:
    """
    Link the position of the candidate PTC and the amino acid.
    """
    aminoacid = []
    for orf, cds, seq in zip(positions_orf, position_cds, seq):
        aminoacid_index = cds // 3
        aminoacid.append({"position": orf, "aminoacid": f"{aminoacid_index+1}{transrate_codon_to_aminoacid(seq)}"})
    return aminoacid


def add_aminoacid_info(candidate_grna: list[dict], aminoacid: list[dict]) -> list[dict]:
    """
    Add the amino acid information to the candidate PTC.
    """
    for cand in candidate_grna:
        for aa in aminoacid:
            if cand["position"] == aa["position"]:
                cand["aminoacid"] = aa["aminoacid"]
    return candidate_grna
