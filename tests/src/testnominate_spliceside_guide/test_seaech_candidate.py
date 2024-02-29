import pytest
from src.nominate_spliceside_guide.search_candidate import search_candidate
from src.create_gene_dataclass import GeneData


def test_no_candidate():
    data = GeneData(
        "NNNNGGNNNNNNNNNNNNNNNNNNNGTNNNCCNNNNNNNNNNNNNNNNNNGANNNNNNNCCNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNN",
        # PAM isn't CC, Donor isn't GT , Acceptor isn't AG
        0,
        100,
        0,
        100,
        3,
        [0, 30, 80],
        [25, 50, 100],
    )
    result = search_candidate(data)
    print(result)
    excepted = ([], [])
    assert result == excepted


def test_has_candidate():
    data = GeneData(
        "NNNNCCNNNNNNNNNNNNNNNNNNNGTNNNCCNNNNNNNNNNNNNNNNNNGTNNNNNNNCCNNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNN",
        0,
        100,
        0,
        100,
        3,
        [0, 30, 80],
        [25, 50, 100],
    )
    result = search_candidate(data)
    print(result)
    excepted = (
        [{"seq": "CCNNNNNNNNNNNNNNNNNAGNN", "exon_num": 3}],
        [
            {"seq": "CCNNNNNNNNNNNNNNNNNNNGT", "exon_num": 1},
            {"seq": "CCNNNNNNNNNNNNNNNNNNGTN", "exon_num": 2},
        ],
    )
    assert result == excepted
