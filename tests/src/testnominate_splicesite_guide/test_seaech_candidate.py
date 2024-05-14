from konezumiaid.nominate_splicesite_guide.search_candidate import search_site_candidate
from konezumiaid.create_gene_dataclass import GeneData


def test_no_candidate():
    data = GeneData(
        "NNNNGGCCNNNNNNNTTTTNNNNNNGTAGNCCNNNNNNNNNNNNNNNNNNGANNNNNNNCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNN",
        # PAM isn't CC,there is TTTT in the sequence, Donor isn't GT , Acceptor isn't AG, and .
        0,
        100,
        0,
        100,
        3,
        [0, 30, 80],
        [25, 50, 100],
    )

    result = search_site_candidate(data)
    excepted = ([], [])
    assert result == excepted


def test_has_candidate():
    data = GeneData(
        "NNNNCCCNNNNNNNNNNNNNNNNNNGTNNNCCNNNNNNNNNNNNNNNNNNGTNNNNNNNCCNNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNN",
        0,
        100,
        0,
        100,
        3,
        [0, 30, 80],
        [25, 50, 100],
    )

    result = search_site_candidate(data)
    excepted = (
        [{"seq": "CCNNNNNNNNNNNNNNNNNAGNN", "exon_index": 3}],
        [
            {"seq": "CCCNNNNNNNNNNNNNNNNNNGT", "exon_index": 1},
            {"seq": "CCNNNNNNNNNNNNNNNNNNGTN", "exon_index": 1},
            {"seq": "CCNNNNNNNNNNNNNNNNNNGTN", "exon_index": 2},
        ],
    )
    assert result == excepted
