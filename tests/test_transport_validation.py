from crystal_transport.validation import run_validation, validation_passes


def test_uniform_and_laminate_limits() -> None:
    result = run_validation(shape=(12, 12, 12))
    assert validation_passes(result)

