# Independent reference for the custom interface

`wavy_wall.py` uses a finite open section with arc-length coordinate
s in [0, 2π] µm and signed gap g(s) = cos(2s)/4 µm, extruded by 0.75 µm.
The curve can be embedded in a plane by taking its horizontal coordinate
to be the integral of sqrt(1 - (dg/ds)^2). Its arc-length rate is one.
There are eight interior crossings of g = ±1/8 µm. The law is zero between
those thresholds and linear outside, with unequal compressive and adhesive
slopes, so the integrated force is nonzero.

The two positive lobes together give an integral of
sqrt(3)/4 - π/12 for the positive part of g - 1/8, from the sine
antiderivative. Translation by π/2 exchanges the positive and negative
lobes, so the compressive side has the same integral before weighting.
The exact force is therefore 0.75 (2 - 3) (sqrt(3)/4 - π/12) pN. The
reference uses this closed form and does not invoke quadrature or any
implementation callback. The comparison requires no quadrature warning
and absolute force error no greater than 1e-12 pN. It verifies this
specified integral, not physical validity or performance across arbitrary
inputs.
