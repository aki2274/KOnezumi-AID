from src.konezumiaid.create_gene_dataclass import GeneData
from src.konezumiaid.get_exon_range import get_exon_range


def test_get_exon_range():
    data = GeneData(
        "NNNNGGCCNNNNNNNAAAANNNNNNGTAGNCCNNNNNNNNNNNNNNNNNNGANNNNNNNCCNNNNNNNNNNNNNNNNAGNNNNNNNNNCCNNNNNNNNNNNNNNNNNNAGNNNNNNNNCCNNNNNNNNNNNNNNNNGTNNNNNNNNNNNN",
        0,
        150,
        0,
        130,
        4,
        [0, 30, 80, 110, 146],
        [25, 50, 100, 137, 149],
    )
    result = get_exon_range(data)
    excepted = [
        [0, 24],
        [25, 44],
        [45, 64],
        [65, 91],
        [92, 94],
    ]
    assert result == excepted
