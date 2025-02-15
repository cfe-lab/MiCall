from typing import Iterator

import numpy as np

DUPLICATION_FACTOR = 1


def stable_random_distribution(high: int, seed: int = 42) -> Iterator[int]:
    if high <= 0:
        return

    rng = np.random.default_rng(seed)
    block = np.arange(high)
    population = np.concatenate([block] * DUPLICATION_FACTOR, axis=0)

    assert len(population) == DUPLICATION_FACTOR * len(block)

    while True:
        index = rng.choice(population)
        yield index
        population[index] = rng.integers(low=0, high=high)
