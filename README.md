# Geometry-resolved adhesive contact: code and numerical results

The CSV and JSON files are the saved results reported in the manuscript. Files 42, 45, 48 and 51 generate the inclined-wall, circular quadrature, circular closed-form and noncircular comparisons, respectively; 54 checks the noncircular boundary roots and circular limit. Protocols 41, 44, 47, 50 and 53 describe the corresponding runs. These synthetic configurations are not experimental observations.

Original results used NumPy 2.3.5, SciPy 1.16.3 and mpmath 1.3.0; the figure script used R 4.5, ggplot2 4.0.2, cowplot 1.2.0 and svglite 2.2.2. Output timings depend on machine and workload. A new run will overwrite saved CSV/JSON files: copy this directory before running the scripts. The reference is evaluated in high precision and rounded to binary64 before the absolute errors are calculated.

Suggested order in the copied directory: python3 42_indicador_geometria.py; python3 45_indicador_gk.py; python3 48_forma_fechada.py; python3 51_teste_nao_circular.py; python3 54_verificar_limite_raizes.py. The Python scripts have fixed cases and write to their own directory. Run Rscript figures_r/render_contact_figures.R from this directory to regenerate the figures from the saved CSV files; only the input path has been adapted for this package.

The review_checks subdirectory contains the post hoc direct-reference checks, isolated-breakpoint diagnostics and fixed Gauss–Legendre comparators; their saved CSV and JSON files are provided. These tests were added after inspection of the primary results and are exploratory.

Author: Renato Socodato, University of Porto. This repository is private during manuscript preparation. Code and results are available from the corresponding author upon reasonable request. The author intends to make the repository public after submission; no public release is represented as having occurred. No licence has yet been selected.

The computed quantities describe idealized contact models and synthetic configurations; they are not measurements from biological specimens. The 52 circular and 18 noncircular cases, near-endpoint selection and post hoc checks are described in the manuscript and protocols. The strongest observed accuracy violations occur in strongly interpenetrating penalty-model configurations.

The 14 individual figure panels and four composite figures are in `figures_r/`; the numerical scripts save their output beside the scripts, so use a copy for re-execution. The scripts in `review_checks/` use the saved primary tables and include exploratory checks.
