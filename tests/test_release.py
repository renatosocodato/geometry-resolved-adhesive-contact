"""Compare the installed interface to the numerical evidence shipped with the paper."""

import csv
import math
from pathlib import Path
import unittest

from geometry_resolved_adhesive_contact import contact_force

ROOT = Path(__file__).resolve().parents[1]


class PaperParity(unittest.TestCase):
    def test_archived_contact_cases(self):
        for circular, filename in [(True, "45_indicador_gk_resultados.csv"),
                                   (False, "51_nao_circular_resultados.csv")]:
            with (ROOT / filename).open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 156 if circular else 54)
            for row in rows:
                with self.subTest(id=row["id"], method=row["method"]):
                    geometry = ({"radius_um": float(row["r"])} if circular else
                                {"rx_um": .36, "ry_um": .28209479177387814, "q_um": .045})
                    value = contact_force(height_um=float(row["h"]), method=row["method"], **geometry)
                    self.assertTrue(math.isclose(value.force_pn, float(row["force"]),
                                                 abs_tol=1e-13, rel_tol=1e-12))
                    ref = float(row["true_reference" if circular else "reference"])
                    exceeded = abs(value.force_pn - ref) > 1e-8
                    self.assertEqual(exceeded, row["failed"] == "True")
                    self.assertEqual(value.integrand_evaluations,
                                     int(row["evaluations" if circular else "law_evals"]))
                    self.assertEqual(value.gap_search_evaluations,
                                     0 if circular else int(row["root_geometry_evals"]))
                    self.assertEqual(bool(value.quadrature_warning), row["warning"] == "True")

    def test_outside_contact_and_invalid_geometry(self):
        result = contact_force(height_um=1, radius_um=.36)
        self.assertEqual(result.force_pn, 0)
        self.assertEqual(result.transition_angles_rad, ())
        with self.assertRaises(ValueError):
            contact_force(height_um=.2, rx_um=.36, ry_um=.1, q_um=.05)


if __name__ == "__main__":
    unittest.main()
