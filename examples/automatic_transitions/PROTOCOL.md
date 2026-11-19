# Automatic monotonicity scan

Use 32 uniform scan intervals and 1e-14 parameter tolerance. Scan the
provided analytical gap derivative and refine slope sign changes. Use the
detected cuts to bracket every declared traction-law threshold. Compare
the computed force with the existing independent closed forms at 1e-12 pN.

Reuse the narrow benchmark's three centers, seven positive depths, exact
tangency and positive clearance per center. Add the wavy wall with two
constitutive thresholds, three extrema and eight crossings. Require the
expected number and locations of cuts and transitions, no quadrature
warning, and force error within tolerance. Count actual slope callbacks,
gap calls during the root search, and force integrand evaluations
separately. Their sum is an unweighted call count, not elapsed time.
Count all calls from monotonicity detection through integration.

These are regression cases with known analytical references. The scan
cannot certify arbitrary continuous slopes: multiple extrema can sit
between two nodes of the same sign. Do not infer a universal guarantee
from these examples.
