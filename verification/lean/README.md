# Lean 4 verification

This module accompanies *Geometry-resolved integration of adhesive contact
forces on extruded surfaces*. It formalizes the real-valued contact model
used in the [paper reproduction snapshot][paper-commit]. You can run the
proofs independently of the Python package and R figure scripts.

[paper-commit]: https://github.com/renatosocodato/geometry-resolved-adhesive-contact/tree/e69b5ef06e64640bf736b7dec2173bb3cd92047d

## Check the proofs

With Elan available on your PATH, enter this directory. The
`lean-toolchain` file selects Lean 4.32.0 and `lake-manifest.json` pins
Mathlib and its transitive dependencies. The verification script also
requires Bash and Python 3.

```sh
cd verification/lean
lake exe cache get \
  Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse \
  Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv \
  Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus \
  Mathlib.Analysis.Calculus.ParametricIntegral \
  Mathlib.Analysis.Calculus.LHopital \
  Mathlib.Analysis.Calculus.ContDiff.Deriv \
  Mathlib.Tactic
bash verify.sh
```

The first command after `cd` downloads the compiled dependencies.
The script checks all 76 theorems in [Contact.lean](Contact.lean) and
records the kernel output in `check.log`. It checks the axiom dependencies
of every theorem against `propext`, `Classical.choice` and `Quot.sound`.

## Statements covered

| Paper result | Representative Lean theorems |
| --- | --- |
| Circular gap and threshold locations, Eqs. (1) and (4) | `circular_interior_root`, `circular_unique`, `circular_endpoint_cases` |
| Contact potential and slope changes, Eq. (2) | `potential_C1`, `actual_circular_jump_at_zero`, `actual_circular_jump_at_range` |
| Force as the negative energy derivative and its circular closed form | `translation_energy_derivative`, `circular_force_closed_form`, `noncircular_energy_force` |
| Noncircular metric and threshold crossing, Eqs. (5) and (6) | `noncircular_metric_is_speed`, `noncircular_monotone`, `noncircular_bracketed_root` |
| Inclined-wall gap and area, §2.2 | `inclined_gap_and_area`, `inclined_surface_area_from_map`, `inclined_physical_axial_cuts` |
| Local force correction at a nondegenerate circular endpoint, §3 | `endpoint_width_delta_limit`, `actual_patch_force_limit`, `actual_patch_force_coefficient_nonzero` |

The proofs apply to the prescribed geometries and constitutive law, under
the hypotheses stated in each theorem. The endpoint result is a one-sided
limit for an isolated region and uses the magnitude of the endpoint gap
excursion. Floating-point discrepancies, quadrature performance and
figures are addressed by the separate [numerical reproduction workflow](../../README.md#reproduce-the-paper).

The Lean code and verification script use the repository's [MIT license](../../LICENSE).
