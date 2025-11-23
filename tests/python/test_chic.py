import chic

def test_add():
    ex = chic.Example()
    assert ex.add(2, 3) == 5
