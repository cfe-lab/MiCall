from typing import Iterator, Optional, Sequence

import random
import numpy as np


def stable_random_distribution(high: int,
                               rng: Optional[random.Random] = None,
                               ) -> Iterator[int]:

    if high <= 0:
        return

    if rng is None:
        rng = random.Random()

    population = np.arange(high)
    weights = np.zeros(high) + 16.5

    while True:
        pweights: Sequence[float] = weights  # type: ignore
        index = rng.choices(population, weights=pweights)[0]
        yield index
        weights[index] *= 0.5
        if weights[index] < 1.0:
            weights += 1.0
