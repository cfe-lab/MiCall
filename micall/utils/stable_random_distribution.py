from typing import Iterator

import numpy as np


def stable_random_distribution(maximum: int) -> Iterator[int]:
    n = maximum

    weights = np.zeros(n) + 1
    forward = np.arange(1, n + 1)
    backwards = np.arange(n, 0, -1)

    while True:
        probabilities = weights / weights.sum()
        index = np.random.choice(n, p=probabilities)
        yield index
        weights[:index] += forward[:index]
        weights[index:] += backwards[index:]
