from __future__ import annotations
from src.create_gene_dataclass import GeneData


def search_splice_candidate(
    ds: GeneData,
) -> tuple[list[dict[int, str], list[dict[int, str]]]]:
    acceptor_candidates = [
        {
            "seq": ds.orf_seq[s - 22 : s + 3][
                ds.orf_seq[s - 22 : s + 3]
                .upper()
                .index("CC") : ds.orf_seq[s - 22 : s + 3]
                .upper()
                .index("CC")
                + 23
            ],
            "exon_num": i + 2,
        }
        for i, s in enumerate(ds.exon_start_list[1:])
        if "CC" in ds.orf_seq[s - 22 : s + 3].upper()[:4]
        and "AG" in ds.orf_seq[s - 2 : s].upper()
    ]

    donor_candidates = [
        {
            "seq": ds.orf_seq[e - 21 : e + 4][
                ds.orf_seq[e - 21 : e + 4]
                .upper()
                .index("CC") : ds.orf_seq[e - 21 : e + 4]
                .upper()
                .index("CC")
                + 23
            ],
            "exon_num": i + 1,
        }
        for i, e in enumerate(ds.exon_end_list[:-1])
        if "CC" in ds.orf_seq[e - 21 : e + 4].upper()[:4]
        and "GT" in ds.orf_seq[e : e + 2].upper()
    ]

    return acceptor_candidates, donor_candidates
