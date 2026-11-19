# Narrow contact and transition-search costs

Run with the installed package and the dependencies in requirements.txt.
From the project root:

```sh
python examples/narrow_contact/run.py --output narrow-contact-recomputed.csv
```

The three strategies use the same [analytical reference](PROTOCOL.md)
and the same 1e-12 pN absolute force tolerance. The [full table](results.csv)
contains 81 runs, including exact-tangency and positive-clearance controls.

| Across 21 active cases | No breakpoints | Analytical breakpoints | Numerically located breakpoints |
|---|---:|---:|---:|
| Above tolerance | 10 | 0 | 0 |
| Maximum absolute force error (pN) | 1.6875e-5 | 2.0470e-18 | 2.0470e-18 |
| Maximum quadrature estimate (pN) | 9.3025e-13 | 1.8735e-16 | 1.8735e-16 |
| Gap evaluations during search | 0 | 0 | 827 |
| Force-integrand evaluations | 10,563 | 1,323 | 1,323 |

All ten unsegmented tolerance failures reported zero estimated error and
no warning. The maximum discrepancy from an analytical root was
2.0e-15 µm. In the located strategy, each zero-force control used three
gap evaluations. Tangency returned the known vertex as one breakpoint
and used 42 integrand calls; clearance returned none and used 21.
The old wavy-wall example also passes using numerically located
transitions at both thresholds.

Search counts include root refinement and endpoint evaluation, with
reuse of gap values across brackets and levels. They exclude the cost
of supplying the monotone partition. An integrand call also evaluates
the metric and the traction, so the counts do not measure wall time.
These results cover the specific configurations here, not a general
error bound or a guarantee for unknown extrema.
