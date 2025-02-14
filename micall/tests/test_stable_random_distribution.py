
from micall.utils.stable_random_distribution import stable_random_distribution


def test_indices_in_range():
    """Test that each index generated is within the range [0, maximum)."""

    maximum = 10
    gen = stable_random_distribution(maximum, seed=123)
    # Grab a bunch of values from the infinite generator

    for _ in range(1000):
        idx = next(gen)
        assert 0 <= idx < maximum, f"Index {idx} out of range [0,{maximum})"
