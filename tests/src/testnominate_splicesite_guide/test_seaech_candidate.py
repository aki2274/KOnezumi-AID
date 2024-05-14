from src.konezumiaid.nominate_splicesite_guide.search_candidate import (
    search_site_candidate,
)
from src.konezumiaid.create_gene_dataclass import GeneData


def test_no_candidate():
    data = GeneData(
        "NNNNGGCCNNNNNNNAAAANNNNNNGTAGNCCNNNNNNNNNNNNNNNNNNGANNNNNNNCCNNNNNNNNNNNNNNNNAGNNNNNNNNNCCNNNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        # PAM isn't CC,there is TTTT(AAAA, because this gRNA is revers compliment), Donor isn't GT , Acceptor isn't AG, and the exon length is a multiple of 3
        0,
        140,
        0,
        140,
        4,
        [0, 30, 80, 110],
        [25, 50, 100, 140],
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
