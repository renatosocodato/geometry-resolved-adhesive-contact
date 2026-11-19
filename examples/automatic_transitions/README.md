# Automatically located monotonicity cuts

Run `python examples/automatic_transitions/run.py --output automatic.csv`
after installing the wheel and requirements.txt. The
[protocol](PROTOCOL.md) specifies the scan and analytical references.
The [results](results.csv) include the error estimated by the quadrature
and the discrepancy from the reference. All runs met 1e-12 pN.

| Cases | Slope calls | Gap calls during search | Integrand calls | Total calls | Maximum actual error (pN) |
|---|---:|---:|---:|---:|---:|
| 21 active narrow contacts | 707 | 827 | 1323 | 2857 | 2.0470e-18 |
| Six tangent or clear controls | 202 | 18 | 189 | 409 | 0 |
| Wavy wall, eight transitions | 36 | 69 | 189 | 294 | 1.2218e-16 |

For the active narrow cases, the same quadrature used 1,323 calls
when both transition points were supplied analytically, 2,150 combined
calls when the monotonicity cuts were supplied, and 10,563 calls
without breakpoints, with ten tolerance failures in the last case.
Counts treat calls equally and are not an elapsed-time comparison.

The scan requires a continuous derivative callback. It locates resolved
slope sign changes, but cannot certify that there are no hidden extrema
in one sampling cell. For known extrema or when isolated-root guarantees
matter, use a partition justified from the geometry.
