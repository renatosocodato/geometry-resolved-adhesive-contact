"""A prescribed nonmonotone wall and a custom traction law with exact force."""

from dataclasses import asdict
import json
import math

from geometry_resolved_adhesive_contact import integrate_contact


def gap_um(s):
    return 0.25 * math.cos(2 * s)


def arc_length_rate_um(s):
    return 1.0  # s is the numerical arc length in micrometres


def traction_pn_um2(gap):
    if gap < -0.125:
        return 2 * (-gap - 0.125)
    if gap > 0.125:
        return -3 * (gap - 0.125)
    return 0.0


TRANSITIONS = tuple(k * math.pi / 6 for k in (1, 2, 4, 5, 7, 8, 10, 11))
LENGTH_UM = 0.75
REFERENCE_PN = -LENGTH_UM * (math.sqrt(3) / 4 - math.pi / 12)
ABS_TOL_PN = 1e-12


def solve():
    return integrate_contact(
        gap_um=gap_um, arc_length_rate_um=arc_length_rate_um,
        traction_pn_um2=traction_pn_um2, interval=(0.0, 2 * math.pi),
        length_um=LENGTH_UM, breakpoints=TRANSITIONS,
        tolerance_pn=ABS_TOL_PN,
    )


if __name__ == "__main__":
    result = solve()
    print(json.dumps({**asdict(result), "reference_pn": REFERENCE_PN,
                      "absolute_error_pn": abs(result.force_pn - REFERENCE_PN),
                      "acceptance_pn": ABS_TOL_PN}, indent=2))
