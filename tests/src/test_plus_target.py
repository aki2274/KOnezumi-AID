import pytest
from src import plus_get_target_num

def test_get_CtoT_target():
    test_seq_list = ['CCANNNNNNNNNNNNNNNNNNTGG','CAGNNNNNNNNNNNNNNNNGG','CAGNNNNNNNNNNNNNNNNGGCGANNNNNNNNNNNNNNNNGG',
                       'CAGCAANNNNNNNNNNNNNGGTGG']
    result = []
    for i in range(len(test_seq_list)):
        seq = test_seq_list[i]
        target = plus_get_target_num.get_CtoT_target(seq)
        result.append(plus_get_target_num.CtoT_target_start(target,seq))
    expected = [[],[0],[0,21],[0,3]]
    assert result == expected


def test_get_AtoG_target():#listの0要素目はマッチするべきだが、配列の一番最後である場合なので検索できなくてもOK？
    test_seq_list = ['CCANNNNNNNNNNNNNNTGG','CCANNNNNNNNNNNNNNTGGNNNN','CCANNNNNNNNNNNNNNNNNTGGCCANNNNNNNNNNNNNNNNNTGG',
                       'CCNCCCNNNNNNNNNNNNNNTGGTGGNNNN']
    result = []
    for i in range(len(test_seq_list)):
        seq = test_seq_list[i]
        target = plus_get_target_num.get_AtoG_target(seq)
        result.append(plus_get_target_num.AtoG_target_start(target,seq))
    expected = [[17],[17],[20,43],[20,23]]
    assert result == expected