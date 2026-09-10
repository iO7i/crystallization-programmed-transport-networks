import numpy as np

from crystal_transport.generators import canonical_morphologies, lamellae
from crystal_transport.morphology import Morphology


def test_morphology_digest_is_reproducible() -> None:
    first = lamellae((16, 16, 16), fraction=0.4)
    second = lamellae((16, 16, 16), fraction=0.4)
    assert first.content_digest() == second.content_digest()
    assert np.array_equal(first.phase, second.phase)


def test_canonical_morphologies_match_phase_fraction() -> None:
    morphologies = canonical_morphologies((16, 16, 16), fraction=0.4)
    fractions = [morphology.volume_fraction for morphology in morphologies]
    assert max(fractions) - min(fractions) == 0.0


def test_serialization_round_trip(tmp_path) -> None:
    original = lamellae((8, 8, 8), fraction=0.5)
    path = tmp_path / "morphology.npz"
    original.save(path)
    restored = Morphology.load(path)
    assert restored.content_digest() == original.content_digest()

