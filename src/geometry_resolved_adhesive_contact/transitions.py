"""Locate gap thresholds on sections partitioned into monotone intervals."""

from collections.abc import Callable, Iterable
from dataclasses import dataclass
import math

from scipy.optimize import brentq


@dataclass(frozen=True)
class TransitionResult:
    breakpoints: tuple[float, ...]
    gap_evaluations: int


def locate_transitions(
    *,
    gap_um: Callable[[float], float],
    interval: tuple[float, float],
    monotone_cuts: Iterable[float],
    thresholds_um: Iterable[float],
    parameter_xtol: float = 1e-14,
) -> TransitionResult:
    """Find threshold locations using a caller-supplied monotone partition.

    Provide interior cuts so that the gap is continuous and monotone in
    each resulting interval. The locator does not verify monotonicity or
    discover extrema. It includes interior cut points exactly on a threshold,
    excludes the overall endpoints, and rejects intervals plateauing at a
    threshold. Tolerance is in section-parameter units, not force units.
    The count is the number of actual gap callback invocations, including
    interval endpoints and root refinement, with reuse across thresholds.
    """
    start, end = interval
    if not (math.isfinite(start) and math.isfinite(end) and start < end):
        raise ValueError("interval must be finite and increasing")
    if not math.isfinite(parameter_xtol) or parameter_xtol <= 0:
        raise ValueError("parameter_xtol must be positive and finite")
    cuts = tuple(sorted(set(monotone_cuts)))
    if any(not math.isfinite(cut) or not start < cut < end for cut in cuts):
        raise ValueError("monotone cuts must be finite and strictly inside interval")
    thresholds = tuple(sorted(set(thresholds_um)))
    if any(not math.isfinite(level) for level in thresholds):
        raise ValueError("thresholds must be finite")
    if not thresholds:
        return TransitionResult((), 0)

    samples = {}

    def gap(coordinate):
        if coordinate not in samples:
            value = gap_um(coordinate)
            if not math.isfinite(value):
                raise ValueError("gap must be finite throughout the search")
            samples[coordinate] = value
        return samples[coordinate]

    bounds = (start, *cuts, end)
    roots = set()
    for left, right in zip(bounds, bounds[1:]):
        gl, gr = gap(left), gap(right)
        for level in thresholds:
            if gl == gr == level:
                raise ValueError("a monotone interval lies on a threshold throughout")
            if gl == level and left != start:
                roots.add(left)
            if gr == level and right != end:
                roots.add(right)
            if gl < level < gr or gr < level < gl:
                root = brentq(lambda s: gap(s) - level, left, right, xtol=parameter_xtol)
                if not left < root < right:
                    raise ArithmeticError("a strict interior transition rounded to a bracket endpoint")
                roots.add(root)
    return TransitionResult(tuple(sorted(roots)), len(samples))
