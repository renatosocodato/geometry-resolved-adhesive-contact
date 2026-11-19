"""Validate the automatic locator against the existing analytical references."""

import argparse
import csv
import math
from pathlib import Path
import sys

import mpmath as mp

from geometry_resolved_adhesive_contact import (
    detect_monotone_cuts, integrate_contact, locate_transitions,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from narrow_contact.run import CENTERS, CURVATURE, HALF_WIDTHS, LENGTH, STIFFNESS, TOLERANCE
from wavy_wall import LENGTH_UM, TRANSITIONS, gap_um, traction_pn_um2


def calculate(gap, slope, traction, end, length, thresholds, cuts, roots, reference):
    scan = detect_monotone_cuts(gap_slope=slope, interval=(0, end))
    found = locate_transitions(gap_um=gap, interval=(0, end), monotone_cuts=scan.cuts,
                               thresholds_um=thresholds)
    assert len(scan.cuts) == len(cuts)
    assert all(abs(a-b) <= 1e-14 for a, b in zip(scan.cuts, cuts))
    assert len(found.breakpoints) == len(roots)
    assert all(abs(a-b) <= 1e-14 for a, b in zip(found.breakpoints, roots))
    force = integrate_contact(
        gap_um=gap, arc_length_rate_um=lambda s: 1.0, traction_pn_um2=traction,
        interval=(0, end), length_um=length, breakpoints=found.breakpoints,
        tolerance_pn=TOLERANCE,
    )
    actual = abs(mp.mpf(force.force_pn) - reference)
    assert actual <= TOLERANCE and not force.quadrature_warning
    return dict(
        reference_pn=mp.nstr(reference, 30), force_pn=force.force_pn,
        abs_error_pn=float(actual), estimate_pn=force.estimated_error_pn,
        slope_evaluations=scan.slope_evaluations, gap_search_evaluations=found.gap_evaluations,
        integrand_evaluations=force.integrand_evaluations,
        total_calls=scan.slope_evaluations + found.gap_evaluations + force.integrand_evaluations,
        cuts=len(scan.cuts), breakpoints=len(found.breakpoints),
    )


def main(output):
    rows = []
    with mp.workdps(70):
        for center in CENTERS:
            for depth in [*(CURVATURE*w**2 for w in HALF_WIDTHS), 0.0, -0.000001]:
                width = mp.sqrt(mp.mpf(depth)/CURVATURE) if depth > 0 else mp.mpf(0)
                roots = (float(mp.mpf(center)-width), float(mp.mpf(center)+width)) if width else (
                    (center,) if depth == 0 else ()
                )
                ref = mp.mpf(4)/3 * LENGTH * STIFFNESS * mp.mpf(depth) * width
                metrics = calculate(
                    lambda s: CURVATURE*(s-center)**2-depth,
                    lambda s: 2*CURVATURE*(s-center),
                    lambda g: STIFFNESS*max(-g, 0),
                    1, LENGTH, (0,), (center,), roots, ref,
                )
                rows.append(dict(case="narrow", center_um=center, depth_um=depth, **metrics))
        metrics = calculate(
            gap_um, lambda s: -0.5*math.sin(2*s), traction_pn_um2,
            2*math.pi, LENGTH_UM, (-0.125, 0.125),
            (math.pi/2, math.pi, 3*math.pi/2), TRANSITIONS,
            -mp.mpf(LENGTH_UM)*(mp.sqrt(3)/4-mp.pi/12),
        )
        rows.append(dict(case="wavy", center_um="", depth_um="", **metrics))
    assert len(rows) == 28
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
