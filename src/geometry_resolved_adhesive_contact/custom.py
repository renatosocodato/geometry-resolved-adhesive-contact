"""Integrate a user-defined plane-normal traction over a prescribed section."""

from collections.abc import Callable, Iterable
from dataclasses import dataclass
import math

from scipy.integrate import quad


@dataclass(frozen=True)
class CustomForceResult:
    force_pn: float
    estimated_error_pn: float
    integrand_evaluations: int
    supplied_breakpoints: tuple[float, ...]
    quadrature_warning: str | None


def integrate_contact(
    *,
    gap_um: Callable[[float], float],
    arc_length_rate_um: Callable[[float], float],
    traction_pn_um2: Callable[[float], float],
    interval: tuple[float, float],
    length_um: float,
    breakpoints: Iterable[float] = (),
    tolerance_pn: float = 1e-8,
    max_subintervals: int = 100,
) -> CustomForceResult:
    """Return the signed plane-normal force on an extruded section.

    The interval covers the entire section to be integrated, with no implicit
    symmetry multiplier. The geometric callbacks return signed gap and arc
    length per unit parameter. The law takes signed gap and returns signed
    plane-normal traction, positive away from the plane. Supply all known
    interior branch transitions or other integrand breakpoints in parameter
    coordinates. The quadrature estimate is not a certified error bound.
    """
    start, end = interval
    if not (math.isfinite(start) and math.isfinite(end) and start < end):
        raise ValueError("interval must be finite and increasing")
    if not (math.isfinite(length_um) and length_um > 0
            and math.isfinite(tolerance_pn) and tolerance_pn > 0):
        raise ValueError("length_um and tolerance_pn must be positive and finite")
    integral_tolerance = tolerance_pn / length_um
    if not math.isfinite(integral_tolerance) or integral_tolerance <= 0:
        raise ValueError("tolerance_pn / length_um is outside the representable range")
    if (not isinstance(max_subintervals, int) or isinstance(max_subintervals, bool)
            or max_subintervals < 1):
        raise ValueError("max_subintervals must be a positive integer")

    points = tuple(sorted(set(breakpoints)))
    if any(not math.isfinite(p) or not start < p < end for p in points):
        raise ValueError("breakpoints must be finite and strictly inside interval")
    if len(points) + 1 > max_subintervals:
        raise ValueError("max_subintervals must cover the intervals defined by the breakpoints")

    def integrand(parameter):
        gap = gap_um(parameter)
        metric = arc_length_rate_um(parameter)
        if not math.isfinite(gap) or not math.isfinite(metric) or metric < 0:
            raise ValueError("geometry must return a finite gap and a finite nonnegative arc-length rate")
        traction = traction_pn_um2(gap)
        value = traction * metric
        if not math.isfinite(traction) or not math.isfinite(value):
            raise ValueError("traction and the weighted integrand must be finite")
        return value

    result = quad(integrand, start, end, epsabs=integral_tolerance, epsrel=0.0,
                  points=points or None, limit=max_subintervals, full_output=1)
    value, estimate, info = result[:3]
    force, error = length_um * value, length_um * estimate
    if not math.isfinite(force) or not math.isfinite(error):
        raise ArithmeticError("integrated force or error estimate is not finite")
    return CustomForceResult(force, error, info["neval"], points,
                             result[3] if len(result) > 3 else None)
