"""Check transition search on a multi-lobe profile and at tangency."""

import math
import unittest

from examples.wavy_wall import ABS_TOL_PN, LENGTH_UM, REFERENCE_PN, TRANSITIONS, gap_um, traction_pn_um2
from geometry_resolved_adhesive_contact import integrate_contact, locate_transitions


class LocateTransitions(unittest.TestCase):
    def test_wavy_wall(self):
        calls = []

        def counted_gap(s):
            calls.append(s)
            return gap_um(s)

        result = locate_transitions(
            gap_um=counted_gap, interval=(0, 2 * math.pi),
            monotone_cuts=(math.pi / 2, math.pi, 3 * math.pi / 2),
            thresholds_um=(-0.125, 0.125),
        )
        self.assertEqual(len(result.breakpoints), len(TRANSITIONS))
        for actual, expected in zip(result.breakpoints, TRANSITIONS):
            self.assertAlmostEqual(actual, expected, delta=1e-14)
        self.assertEqual(len(calls), len(set(calls)))
        self.assertEqual(result.gap_evaluations, len(calls))
        force = integrate_contact(
            gap_um=gap_um, arc_length_rate_um=lambda s: 1.0,
            traction_pn_um2=traction_pn_um2, interval=(0, 2 * math.pi),
            length_um=LENGTH_UM, breakpoints=result.breakpoints,
            tolerance_pn=ABS_TOL_PN,
        )
        self.assertAlmostEqual(force.force_pn, REFERENCE_PN, delta=ABS_TOL_PN)
        self.assertIsNone(force.quadrature_warning)

    def test_tangency_and_global_endpoints(self):
        result = locate_transitions(
            gap_um=lambda s: s*s, interval=(-1, 1), monotone_cuts=(0,),
            thresholds_um=(0, 1),
        )
        self.assertEqual(result.breakpoints, (0,))
        self.assertEqual(result.gap_evaluations, 3)

    def test_plateau_at_threshold(self):
        with self.assertRaisesRegex(ValueError, "lies on a threshold"):
            locate_transitions(gap_um=lambda s: 0, interval=(0, 1),
                               monotone_cuts=(), thresholds_um=(0,))
