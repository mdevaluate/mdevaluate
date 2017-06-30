from pytest import approx

from mdevaluate import pbc
import numpy as np


def test_pbc_diff():
    x = np.random.rand(10, 3)
    y = np.random.rand(10, 3)
    box = np.ones((3,))

    assert (pbc.pbc_diff(x, x, box) == approx(0)).all()
    dxy = (pbc.pbc_diff(x, y, box)**2).sum(axis=1)**0.5
    assert (dxy <= 0.75**0.5).all()
