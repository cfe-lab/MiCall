
from micall.utils.stable_random_distribution import stable_random_distribution
import numpy as np
from itertools import islice


def test_indices_in_range():
    """Test that each index generated is within the range [0, maximum)."""

    maximum = 10
    gen = stable_random_distribution(maximum, seed=123)
    # Grab a bunch of values from the infinite generator

    for _ in range(1000):
        idx = next(gen)
        assert 0 <= idx < maximum, f"Index {idx} out of range [0,{maximum})"


def test_bounds_are_reachable():
    """Test that both min and max-1 can be generated."""

    maximum = 999
    gen = stable_random_distribution(maximum, seed=123)
    lst = islice(gen, 1000)

    assert 0 in lst
    assert (maximum-1) in lst


def test_deterministic_output_with_seed():
    """
    Test that the generator produces the same sequence when
    re-seeded with the same seed.
    """

    maximum = 15
    seed = 456
    gen1 = stable_random_distribution(maximum, seed=seed)
    gen2 = stable_random_distribution(maximum, seed=seed)

    # Compare the first 50 generated values.
    values1 = [next(gen1) for _ in range(50)]
    values2 = [next(gen2) for _ in range(50)]
    assert values1 == values2, \
        "Generators with the same seed produced different outputs."


def test_different_seeds_differ():
    """
    A sanity check that different seeds usually lead to a different sequence.
    """

    maximum = 15
    gen1 = stable_random_distribution(maximum, seed=789)
    gen2 = stable_random_distribution(maximum, seed=987)

    # Compare the first 50 generated values: while not guaranteed to
    # be different, it is extremely unlikely that the two sequences
    # are identical.
    values1 = [next(gen1) for _ in range(50)]
    values2 = [next(gen2) for _ in range(50)]

    assert values1 != values2, \
        "Generators with different seeds produced identical sequences."


def test_fair_distribution_behavior():
    """
    Test that the stable_random_distribution leads to outputs that are
    more 'spread out' than a simple uniform generator.

    Idea:
      - Generate a long sequence from our generator.
      - Compute the average absolute difference (jump) between indices.
      - Do the same for a uniformly random generator over the same range.
      - With the adaptive update, values should tend to be farther apart.
    """

    maximum = 1_000
    num_samples = 3_000
    for seed in range(100):
        # Gather samples from our generator.
        gen = stable_random_distribution(maximum, seed=seed)
        samples = np.array([next(gen) for _ in range(num_samples)])
        diff_stable = np.abs(np.diff(np.sort(samples))) ** 2
        avg_diff_stable = diff_stable.mean()

        # For comparison, generate num_samples indices uniformly at random.
        rng = np.random.default_rng(seed)
        uniform_samples = rng.choice(maximum, size=num_samples)
        diff_uniform = np.abs(np.diff(np.sort(uniform_samples))) ** 2
        avg_diff_uniform = diff_uniform.mean()

        # Our expectation: the stable generator should have larger jumps
        # on average. We include a tolerance, because both sequences are
        # random.
        assert avg_diff_stable >= avg_diff_uniform, (
            f"Expected stable generator to have a higher average jump than a uniform generator: "
            f"stable {avg_diff_stable} vs uniform {avg_diff_uniform}"
        )
