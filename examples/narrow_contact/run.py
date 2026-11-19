"""Reproduce the narrow-contact comparison using the installed v0.2.0 wheel."""

import argparse
import csv
from importlib.metadata import version
import math
from pathlib import Path

import mpmath as mp

from geometry_resolved_adhesive_contact import integrate_contact, locate_transitions


CURVATURE = 0.25
LENGTH = 0.75
STIFFNESS = 20.0
TOLERANCE = 1e-12
CENTERS = (0.371234, 0.5, 0.81237)
HALF_WIDTHS = (0.15, 0.05, 0.015, 0.005, 0.001, 0.0002, 0.00005)


def run(output):
    assert version("geometry-resolved-adhesive-contact") == "0.2.0"
    rows = []
    with mp.workdps(70):
        tol = mp.mpf(TOLERANCE)
        for center in CENTERS:
            for depth in [*(CURVATURE * w**2 for w in HALF_WIDTHS), 0.0, -0.000001]:
                width = mp.sqrt(mp.mpf(depth) / CURVATURE) if depth > 0 else mp.mpf(0)
                roots = (float(mp.mpf(center) - width), float(mp.mpf(center) + width)) if width else ()
                ref = mp.mpf(4) / 3 * LENGTH * STIFFNESS * mp.mpf(depth) * width if width else mp.mpf(0)
                def gap(s):
                    return CURVATURE * (s - center)**2 - depth

                located = locate_transitions(
                    gap_um=gap, interval=(0, 1), monotone_cuts=(center,), thresholds_um=(0,),
                )
                expected = roots if width else (center,) if depth == 0 else ()
                assert len(located.breakpoints) == len(expected)
                root_error = max((abs(a - b) for a, b in zip(located.breakpoints, expected)), default=0.0)
                assert root_error <= 1e-14
                strategies = (
                    ("adaptive", (), 0, None),
                    ("segmented", roots, 0, None),
                    ("located", located.breakpoints, located.gap_evaluations, root_error),
                )
                for strategy, points, search_count, root_error_um in strategies:
                    result = integrate_contact(
                        gap_um=gap,
                        arc_length_rate_um=lambda s: 1.0,
                        traction_pn_um2=lambda g: STIFFNESS * max(-g, 0.0),
                        interval=(0, 1), length_um=LENGTH, breakpoints=points,
                        tolerance_pn=TOLERANCE,
                    )
                    actual = abs(mp.mpf(result.force_pn) - ref)
                    estimate = result.estimated_error_pn
                    assert math.isfinite(estimate) and estimate >= 0
                    warning = result.quadrature_warning or ""
                    if strategy == "located":
                        assert actual <= tol and not warning
                    rows.append(dict(
                        center_um=center, depth_um=depth, half_width_um=float(width), strategy=strategy,
                        reference_pn=mp.nstr(ref, 30), force_pn=result.force_pn,
                        abs_error_pn=float(actual), estimate_pn=estimate,
                        search_gap_evaluations=search_count,
                        integrand_evaluations=result.integrand_evaluations, warning=warning,
                        breakpoints=len(points), max_root_error_um=root_error_um,
                        zero_on_active_support=bool(depth > 0 and result.force_pn == 0),
                        above_tolerance=bool(actual > tol),
                        unflagged_tolerance_failure=bool(actual > tol and estimate <= TOLERANCE and not warning),
                    ))
    assert len(rows) == 81
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
