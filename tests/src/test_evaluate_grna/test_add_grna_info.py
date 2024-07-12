from konezumiaid.evaluate_grna.add_grna_info import link_to_crisperdirect


def test_link_to_crisperdirect():
    candidates = [
        {"seq": "ATCGATCGATCGATCGATCGATCG", "exon_index": 1},
        {"seq": "GCTAGCTAGCTAGCTAGCTAGCTA", "exon_index": 2},
    ]
    expected = [
        {
            "seq": "ATCGATCGATCGATCGATCGATCG",
            "exon_index": 1,
            "link_to_crisprdirect": "https://crispr.dbcls.jp/?userseq=ATCGATCGATCGATCGATCGATCG&pam=NGG&db=mm39",
        },
        {
            "seq": "GCTAGCTAGCTAGCTAGCTAGCTA",
            "exon_index": 2,
            "link_to_crisprdirect": "https://crispr.dbcls.jp/?userseq=GCTAGCTAGCTAGCTAGCTAGCTA&pam=NGG&db=mm39",
        },
    ]
    assert link_to_crisperdirect(candidates) == expected
