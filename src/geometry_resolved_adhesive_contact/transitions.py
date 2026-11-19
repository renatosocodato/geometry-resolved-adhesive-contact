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


@dataclass(frozen=True)
class MonotonicityResult:
    cuts: tuple[float, ...]
    slope_evaluations: int


def detect_monotone_cuts(
    *,
    gap_slope: Callable[[float], float],
    interval: tuple[float, float],
    scan_intervals: int = 32,
    parameter_xtol: float = 1e-14,
) -> MonotonicityResult:
    """Find candidate extrema by scanning the supplied continuous gap slope.

    Refine slope sign changes by bracketing, and retain isolated sampled
    zeros. This is a resolution-dependent search, not proof of monotonicity:
    two extrema within one scan cell may go undetected. Consecutive zero
    samples are ambiguous and require explicit geometry information.
    """
    start, end = interval
    if not (math.isfinite(start) and math.isfinite(end) and start < end
            and math.isfinite(end - start)):
        raise ValueError("interval must have finite increasing bounds and width")
    if type(scan_intervals) is not int or scan_intervals < 2:
        raise ValueError("scan_intervals must be an integer of at least two")
    if not (math.isfinite(parameter_xtol) and parameter_xtol > 0):
        raise ValueError("parameter_xtol must be positive and finite")
    nodes = [start + (end - start) * (i / scan_intervals)
             for i in range(scan_intervals)] + [end]
    if any(a >= b for a, b in zip(nodes, nodes[1:])):
        raise ValueError("scan nodes are indistinguishable at this parameter scale")
    samples = {}

    def slope(t):
        if t not in samples:
            value = gap_slope(t)
            if not math.isfinite(value):
                raise ValueError("gap slope must be finite throughout the search")
            samples[t] = value
        return samples[t]

    values = [slope(t) for t in nodes]
    cuts = {t for t, value in zip(nodes[1:-1], values[1:-1]) if value == 0}
    for left, right, sl, sr in zip(nodes, nodes[1:], values, values[1:]):
        if sl == sr == 0:
            raise ValueError("consecutive zero slopes: monotonicity scan is ambiguous")
        if sl < 0 < sr or sr < 0 < sl:
            root = brentq(slope, left, right, xtol=parameter_xtol)
            # A stationary point within tolerance of an internal scan node
            # can round to that node. Only the global endpoints are excluded.
            if start < root < end:
                cuts.add(root)
    return MonotonicityResult(tuple(sorted(cuts)), len(samples))
