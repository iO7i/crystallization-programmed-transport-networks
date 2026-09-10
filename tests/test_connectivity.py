import numpy as np

from crystal_transport.connectivity import analyze_connectivity
from crystal_transport.morphology import Morphology


def test_directional_periodic_connectivity() -> None:
    phase = np.zeros((8, 8, 8), dtype=bool)
    phase[:, 2, 2] = True
    metrics = analyze_connectivity(Morphology(phase, name="x_channel"))
    assert metrics.percolates == (True, False, False)


def test_disconnected_fixture_does_not_percolate() -> None:
    phase = np.zeros((8, 8, 8), dtype=bool)
    phase[2:4, 2:4, 2:4] = True
    metrics = analyze_connectivity(Morphology(phase, name="island"))
    assert not any(metrics.percolates)
    assert metrics.components == 1

