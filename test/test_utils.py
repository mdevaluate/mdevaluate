from copy import copy
import pytest
import numpy as np

from mdevaluate import utils


@pytest.fixture
def logdata(request):
    xdata = np.logspace(-1, 3, 50)
    ydata = np.exp(- (xdata)**0.7)
    return xdata, ydata


def test_filon_fourier_transformation(logdata):
    xdata, ydata = logdata

    xdata_zero = copy(xdata)
    xdata_zero[0] = 0
    _, filon = utils.filon_fourier_transformation(xdata_zero, ydata)
    assert not np.isnan(filon).any(), 'There are NaN values in the filon result!'

    freqs = np.logspace(-4, 1)
    filon_freqs, filon_imag = utils.filon_fourier_transformation(
        xdata, xdata, frequencies=freqs, derivative='linear', imag=True
        )

    assert (freqs == filon_freqs).all()

    freqs, filon_real = utils.filon_fourier_transformation(
        xdata, xdata, frequencies=freqs, derivative='linear', imag=False
        )
    assert np.isclose(filon_imag.real, filon_real).all()


def test_histogram():
    data = np.random.rand(100)
    bins = np.linspace(0, 1)
    np_hist = np.histogram(data, bins=bins)[0]
    ut_hist = utils.histogram(data, bins=bins)[0]
    assert (np_hist == ut_hist).all()

    bins = np.linspace(0.3, 1.5)
    np_hist = np.histogram(data, bins=bins)[0]
    ut_hist = utils.histogram(data, bins=bins)[0]
    assert (np_hist == ut_hist).all()
