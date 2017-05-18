from mdevaluate import atoms


def test_compare_regex():
    assert atoms.compare_regex(['OW', ], 'O')[0] == False
    assert atoms.compare_regex(['WO', ], 'O')[0] == False
    assert atoms.compare_regex(['O', ], 'O')[0] == True
