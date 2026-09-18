"""
Deterministic random-number source.

CarMaker seeds all stochastic sensor effects from a single Test Run parameter
(``RandomSeed``, GUI -> Parameters -> Additional Parameters; see Reference
Manual -> Sensors -> Radar Cross Section, "Note").  We mirror that behaviour:
one seed per sensor instance, reproducible across runs.

Only three primitives are needed by the radar models, all of which have direct
C++ equivalents (``std::mt19937_64`` + ``std::normal_distribution`` etc.), so
the algorithm ports without change.
"""

from __future__ import annotations

import math
import random
from typing import Optional


class Rng:
    """Thin wrapper around a Mersenne-Twister stream."""

    __slots__ = ("_r", "seed")

    def __init__(self, seed: Optional[int] = None):
        self.seed = 0 if seed is None else int(seed)
        self._r = random.Random(self.seed)

    def reset(self, seed: Optional[int] = None) -> None:
        if seed is not None:
            self.seed = int(seed)
        self._r.seed(self.seed)

    # ------------------------------------------------------------ primitives
    def uniform(self, lo: float = 0.0, hi: float = 1.0) -> float:
        return self._r.uniform(lo, hi)

    def normal(self, mean: float = 0.0, sigma: float = 1.0) -> float:
        return self._r.gauss(mean, sigma)

    def exponential(self, mean: float) -> float:
        """Exponentially distributed value with the given mean.

        Used for the Swerling-1 RCS fluctuation (Reference Manual Eq. 509):
        p(sigma) = (1/sigma_bar) * exp(-sigma / sigma_bar).
        """
        if mean <= 0.0:
            return 0.0
        u = self._r.random()
        if u <= 0.0:
            u = 1e-300
        return -mean * math.log(u)

    def poisson(self, lam: float) -> int:
        """Knuth's algorithm - small lambda only (ClutterObjMean is O(1..100))."""
        if lam <= 0.0:
            return 0
        if lam < 30.0:
            limit = math.exp(-lam)
            k = 0
            p = 1.0
            while True:
                p *= self._r.random()
                if p <= limit:
                    return k
                k += 1
        # normal approximation for large lambda
        return max(0, int(round(self._r.gauss(lam, math.sqrt(lam)))))
