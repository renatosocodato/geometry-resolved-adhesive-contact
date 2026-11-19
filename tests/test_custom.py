"""Verify the custom interface against an analytical reference."""

import math
import unittest

from examples.wavy_wall import ABS_TOL_PN, REFERENCE_PN, TRANSITIONS, solve
from geometry_resolved_adhesive_contact import integrate_contact


class CustomContact(unittest.TestCase):
    def test_nonmonotone_wall_against_closed_form(self):
        result = solve()
        self.assertTrue(math.isclose(result.force_pn, REFERENCE_PN,
                                     rel_tol=0, abs_tol=ABS_TOL_PN))
        self.assertIsNone(result.quadrature_warning)
        self.assertEqual(result.supplied_breakpoints, TRANSITIONS)

    def test_section_metric_sign_and_extrusion(self):
        result = integrate_contact(
            gap_um=lambda s: 1.0, arc_length_rate_um=lambda s: s,
            traction_pn_um2=lambda g: -2.0, interval=(0, 3), length_um=2,
        )
        self.assertAlmostEqual(result.force_pn, -18.0)

    def test_endpoint_breakpoint_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "strictly inside"):
            integrate_contact(
                gap_um=lambda s: s, arc_length_rate_um=lambda s: 1.0,
                traction_pn_um2=lambda g: g, interval=(0, 1), length_um=1,
                breakpoints=(0,),
            )
