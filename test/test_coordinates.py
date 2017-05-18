import os
import pytest

import mdevaluate
from mdevaluate import coordinates


@pytest.fixture
def trajectory(request):
    return mdevaluate.open(os.path.join(os.path.dirname(__file__), 'data/water'))


def test_coordinates_getitem(trajectory):
    """
    Tests for the Coordinates class.
    """
    assert isinstance(trajectory[0], coordinates.CoordinateFrame)
    assert isinstance(trajectory[len(trajectory) - 1], coordinates.CoordinateFrame)
    i = 0
    dt = trajectory[1].time - trajectory[0].time
    for f in trajectory:
        assert f.step == i
        assert round(f.time, 3) == round(i * dt, 3)
        i += 1
    sl = trajectory[0::10]
    assert isinstance(sl, coordinates.Coordinates)
    i = 0
    for f in sl:
        assert f.step == i
        assert round(f.time, 3) == round(i * dt, 3)
        i += 10
