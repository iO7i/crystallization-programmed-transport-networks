import numpy as np

from crystal_transport.crystalline_perturbation import (
    baseline_flux_score,
    build_diffusivity_field,
    select_crystalline_like,
)
from crystal_transport.generators import gyroid_level_set
from crystal_transport.transport import (
    solve_effective_diffusivity,
    solve_effective_diffusivity_field,
)


def test_placement_is_exact_and_seeded() -> None:
    morphology = gyroid_level_set((8, 8, 8), fraction=0.4)
    baseline = tuple(
        solve_effective_diffusivity(morphology, axis=axis, d_transport=1.0, d_matrix=1.0e-3)
        for axis in range(3)
    )
    score = baseline_flux_score(morphology, baseline, 1.0, 1.0e-3)
    first = select_crystalline_like(morphology, 0.2, "baseline-flux-ranked", 7, score)
    second = select_crystalline_like(morphology, 0.2, "baseline-flux-ranked", 7, score)
    assert first.selected_voxel_count == round(0.2 * int(np.count_nonzero(morphology.phase)))
    assert first.selection_digest == second.selection_digest
    assert np.all(~first.crystalline_mask | morphology.phase)


def test_equal_mobility_label_is_a_noop() -> None:
    morphology = gyroid_level_set((8, 8, 8), fraction=0.4)
    mask = select_crystalline_like(morphology, 0.2, "random", 7).crystalline_mask
    field = build_diffusivity_field(morphology, mask, 1.0, 1.0, 1.0e-3)
    expected = np.where(morphology.phase, 1.0, 1.0e-3)
    assert np.array_equal(field, expected)
    scalar = solve_effective_diffusivity(morphology, axis=2, d_transport=1.0, d_matrix=1.0e-3)
    voxel = solve_effective_diffusivity_field(morphology, field, axis=2)
    assert np.isclose(voxel.effective_diffusivity, scalar.effective_diffusivity, rtol=1e-10)
