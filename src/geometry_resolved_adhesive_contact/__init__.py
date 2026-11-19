"""Quadrature for the adhesive contact models defined in the companion paper."""

from dataclasses import dataclass
import math

from scipy.integrate import quad
from scipy.optimize import brentq

from .custom import CustomForceResult, integrate_contact

__all__ = ["CustomForceResult", "ForceResult", "contact_force", "integrate_contact"]


@dataclass(frozen=True)
class ForceResult:
    force_pn: float
    estimated_error_pn: float
    integrand_evaluations: int
    gap_search_evaluations: int
    transition_angles_rad: tuple[float, ...]
    supplied_breakpoints_rad: tuple[float, ...]
    quadrature_warning: str | None


def contact_force(*, height_um: float, radius_um: float | None = None,
                  rx_um: float | None = None, ry_um: float | None = None,
                  q_um: float = 0.0, length_um: float = 0.75,
                  stiffness_pn_um3: float = 20.0, adhesion_pn_um: float = 0.002,
                  range_um: float = 0.1, tolerance_pn: float = 1e-8,
                  method: str = "all") -> ForceResult:
    """Integrate the force conjugate to the plane-normal height.

    Supply either a circular radius or both noncircular coefficients rx and ry.
    The noncircular profile uses x=rx*sin(t), g=h+ry*cos(t)+q*cos(2*t),
    with ry > 4*abs(q) to ensure a monotone gap over (0, pi).
    Methods: gk (no supplied roots), guarded (supply all roots if one lies
    within pi/16 of an endpoint), all (supply all interior roots).
    The quadrature estimate is not a certified error bound.
    """
    if method not in {"gk", "guarded", "all"}:
        raise ValueError("method must be gk, guarded, or all")
    circle = radius_um is not None
    if circle:
        if rx_um is not None or ry_um is not None or q_um != 0:
            raise ValueError("give a circular radius or noncircular coefficients")
        dimensions = (radius_um,)
    else:
        if rx_um is None or ry_um is None or ry_um <= 4 * abs(q_um):
            raise ValueError("the noncircular wall requires rx and ry > 4*abs(q)")
        dimensions = (rx_um, ry_um)
    positive = (*dimensions, length_um, stiffness_pn_um3, range_um, tolerance_pn)
    if not all(math.isfinite(x) and x > 0 for x in positive):
        raise ValueError("dimensions, stiffness, range and tolerance must be positive and finite")
    if not all(math.isfinite(x) for x in (height_um, adhesion_pn_um, q_um)) or adhesion_pn_um < 0:
        raise ValueError("height, adhesion and q must be finite, and adhesion nonnegative")

    def gap(theta):
        if circle:
            return height_um + radius_um * math.cos(theta)
        return height_um + ry_um * math.cos(theta) + q_um * math.cos(2 * theta)

    search_evaluations = 0

    def residual(theta, boundary):
        nonlocal search_evaluations
        search_evaluations += 1
        return gap(theta) - boundary

    roots = []
    if method != "gk":
        for boundary in (0.0, range_um):
            if circle:
                cosine = (boundary - height_um) / radius_um
                if -1 < cosine < 1:
                    roots.append(math.acos(cosine))
            elif residual(0.0, boundary) * residual(math.pi, boundary) < 0:
                roots.append(brentq(residual, 0.0, math.pi, args=(boundary,),
                                    xtol=1e-13, rtol=1e-14))
    roots.sort()
    split = method == "all" or (
        method == "guarded" and any(min(t, math.pi - t) < math.pi / 16 for t in roots)
    )
    points = roots if split else None

    def integrand(theta):
        g = gap(theta)
        if g <= 0:
            value = stiffness_pn_um3 * g
        elif g >= range_um:
            value = 0.0
        else:
            value = 4 * adhesion_pn_um * g / range_um**2 * (1 - (g / range_um)**2)
        if circle:
            return value
        arc = math.hypot(rx_um * math.cos(theta),
                         -ry_um * math.sin(theta) - 2 * q_um * math.sin(2 * theta))
        return value * arc

    factor = 2 * radius_um * length_um if circle else 2 * length_um
    result = quad(integrand, 0.0, math.pi, epsabs=tolerance_pn / factor,
                  epsrel=0.0, points=points, limit=100, full_output=1)
    value, estimate, info = result[:3]
    return ForceResult(-factor * value, factor * estimate, info["neval"],
                       search_evaluations, tuple(roots), tuple(points or ()),
                       result[3] if len(result) > 3 else None)
