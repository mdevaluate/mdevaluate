from mdevaluate import checksum

import numpy as np


def test_checksum():
    salt = checksum.SALT
    checksum.SALT = b''
    assert checksum.checksum(1) == 253040863595743782930149460018847280986
    assert checksum.checksum('42') == checksum.checksum(42)
    cs1 = checksum.checksum(999)
    checksum.SALT = b'999'
    assert cs1 != checksum.checksum(999)

    a = np.array([1, 2, 3])
    assert checksum.checksum(a) == checksum.checksum(a.tobytes())

    checksum.SALT = salt


def test_version():

    @checksum.version(1)
    def f1():
        pass

    cs1 = checksum.checksum(f1)

    @checksum.version(1)
    def f1(x, y):
        return x + y

    assert cs1 == checksum.checksum(f1)

    @checksum.version(2)
    def f1(x, y):
        pass

    assert cs1 != checksum.checksum(f1)
