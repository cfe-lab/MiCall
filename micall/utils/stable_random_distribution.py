from typing import Iterator, Optional

import random
import numpy as np

DUPLICATION_FACTOR = 1


def stable_random_distribution(high: int,
                               rng: Optional[random.Random] = None,
                               ) -> Iterator[int]:

    if high <= 0:
        return

    if rng is None:
        rng = random.Random()

    maximum = high - 1
    block = np.arange(high)
    population = np.concatenate([block] * DUPLICATION_FACTOR, axis=0)

    assert len(population) == DUPLICATION_FACTOR * len(block)

    while True:
        choice = rng.randint(0, maximum)
        index = population[choice]
        yield index
        population[index] = rng.randint(0, maximum)
