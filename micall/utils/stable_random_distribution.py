from typing import Iterator, Sequence

import numpy as np
import random


def stable_random_distribution(maximum: int, seed: int = 42) -> Iterator[int]:
    if maximum <= 0:
        return

    n = maximum
    rng = random.Random(seed)

    population = np.arange(n)
    forward = np.arange(1, n + 1)
    backwards = np.copy(np.flip(forward))
    np_weights = np.zeros(n) + 1

    while True:
        top = np.max(np_weights) + 1
        weights: Sequence[float] = top - np_weights  # type: ignore
        index = rng.choices(population=population, weights=weights)[0]
        yield index

        if index == 0:
            np_weights += backwards
        else:
            np_weights[:(index + 1)] += forward[-(index + 1):]
            np_weights[(index + 1):] += backwards[1:-index]

        # Prevent overflow.
        np_weights = np_weights - np_weights.min()
