# Narrow contact and approach to tangency

Test the installed v0.2.0 custom interface on a finite extruded section
parametrized by arc length s in [0,1] µm. Let g(s) = c(s-s₀)² − δ,
c = 1/4 µm⁻¹, extrusion L = 3/4 µm and traction T(g) = K max(-g,0),
K = 20 pN/µm³. This is a unilateral linear penalty test law with adhesion
set to zero to isolate the compact support. The curve is realizable:
|g'| < 1 throughout, so the horizontal derivative can be sqrt(1 − g'²).
The metric is unity. Decreasing δ to zero makes the two transitions
coalesce into tangency.

For δ > 0, w = sqrt(δ/c), the support is (s₀−w,s₀+w), and for all selected
states it lies inside [0,1]. Direct integration of the quadratic traction
gives F = (4/3) L K δ sqrt(δ/c). For δ ≤ 0 the exact force is zero.
Evaluate the same expression at 70 decimal digits using the exact binary64
input coefficients for both quadrature strategies. The reference invokes
neither the numerical integrand nor quadrature.

Fixed matrix: centers 0.371234, 0.5, and 0.81237 µm, seven positive target
half-widths 0.15, 0.05, 0.015, 0.005, 0.001, 0.0002 and 0.00005 µm,
and exact tangency and positive clearance (δ = −0.000001 µm) controls.
The central position tests a case whose minimum is sampled at the interval
midpoint. All are constructed test cases, not a distribution of contacts.

Compare default adaptive integration without breakpoints against the same
integration with the two exact geometric transition locations rounded to
binary64. Both use absolute force tolerance 1e-12 pN and the default
100-subinterval budget. Record measured absolute error, reported estimate,
evaluation count and warning. Count unflagged tolerance failures when the
measured error is above tolerance, the estimate is at or below it and
there is no warning. Record zero returns on positive support separately
even if their absolute error is below tolerance. A pointwise underestimate
alone is descriptive and is not a tolerance failure. Stop after 54 runs.

## Extension: bracketed transition search

For the same cases add a third strategy. Locate g = 0 over two monotone
intervals [0,s₀] and [s₀,1] using the supplied vertex as the cut. Use
1e-14 absolute parameter tolerance for the Brent solves. Count actual
gap calls during search separately from traction integrand evaluations.
For penetration, compare the two numerical roots with the analytical
locations. For exact tangency expect one interior threshold location and
for positive clearance expect none. Require the same force tolerance and
no quadrature warning for the located strategy. This benchmark includes
no cost of inferring the monotone partition, because its vertex is known.
