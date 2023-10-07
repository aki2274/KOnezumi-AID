import pytest
from src import plus_target

def test_get_CtoT_target():
    test_seq_list = ['CCANNNNNNNNNNNNNNNNNNTGG','CAGNNNNNNNNNNNNNNNNGG','CAGNNNNNNNNNNNNNNNNGGCGANNNNNNNNNNNNNNNNGG',
                       'CAGCAANNNNNNNNNNNNNGGTGG']
    result = []
    for i in range(len(test_seq_list)):
        seq = test_seq_list[i]
        target = plus_target.get_CtoT_target(seq)
        result.append(plus_target.CtoT_target_start(target,seq))
    expected = [[],[0],[0,21],[0,3]]
    assert result == expected


def test_get_AtoG_target():#listの0要素目はマッチするべきだが、配列の一番最後である場合なので検索できなくてもOK？
    test_seq_list = ['CCANNNNNNNNNNNNNNTGG','CCANNNNNNNNNNNNNNTGGNNNN','CCANNNNNNNNNNNNNNNNNTGGCCANNNNNNNNNNNNNNNNNTGG',
                       'CCNCCCNNNNNNNNNNNNNNTGGTGGNNNN']
    result = []
    for i in range(len(test_seq_list)):
        seq = test_seq_list[i]
        target = plus_target.get_AtoG_target(seq)
        result.append(plus_target.AtoG_target_start(target,seq))
    expected = [[],[17],[20,43],[20,23]]
    assert result == expected