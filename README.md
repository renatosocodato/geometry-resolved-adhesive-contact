# Geometry-Resolved Adhesive Contact

Version 0.2.0 provides a Python interface for the plane-normal force on
prescribed circular and monotone-gap noncircular extruded walls against a rigid
plane, using the piecewise-smooth adhesive law of the companion paper.
All inputs use micrometres and piconewtons. Stiffness is in pN/µm³ and adhesion
energy per area is in pN/µm. The quantity returned by the quadrature is an error
estimate within the assumed model.

## Install and run

Install from the source directory with `python -m pip install .`, or install the
wheel supplied with the release. Use Python 3.11 or later.

```sh
adhesive-contact --height-um 0.32 --radius-um 0.28209479177387814
adhesive-contact --height-um 0.30 --rx-um 0.36 --ry-um 0.28209479177387814 --q-um 0.045
```

```python
from geometry_resolved_adhesive_contact import contact_force

result = contact_force(height_um=0.32, radius_um=0.28209479177387814)
print(result.force_pn, result.estimated_error_pn)
```

Use `--method gk`, `guarded` or `all` to choose how interior transitions are
supplied. The JSON and Python outputs distinguish force-integrand evaluations
from gap-search evaluations, report the supplied breakpoints, and expose
quadrature warnings. The noncircular model requires `ry > 4*abs(q)`.
The model assumes prescribed walls and does not solve for deformation.

## Supply your own section and traction law

Use the Python API when your geometry or law differs from the built-in
models. Pass the signed gap and arc length per unit parameter as functions
of the section coordinate, and pass the signed plane-normal traction as a
function of gap. Supply any known interior breakpoints in the section
coordinate, including all crossings of thresholds in a nonmonotone gap.
The section integral is multiplied by extrusion length once. Pass the full
section unless you intend to compute the force on only part of it.

```python
from geometry_resolved_adhesive_contact import integrate_contact

result = integrate_contact(
    gap_um=lambda t: t - 0.5,
    arc_length_rate_um=lambda t: 1.0,
    traction_pn_um2=lambda g: 10.0 * max(-g, 0.0),
    interval=(0.0, 1.0),
    length_um=0.75,
    breakpoints=(0.5,),
)
print(result.force_pn)  # 0.9375 pN, positive away from the plane
```

The callback interface leaves transition location to the caller. Check
quadrature warnings and use an independent reference where accuracy matters.
See [the wavy-wall example](examples/README.md) for a nonmonotone geometry
with eight transitions and a closed-form force. Run it from this archive
with `python examples/wavy_wall.py` after installing the package.

Version 0.1.0 is the original software snapshot for the paper. The original
benchmark scripts and archived evidence are also included here unchanged.
The CLI continues to accept the built-in geometries and contact law.

## Reproduce the paper

Run from the unpacked **paper reproduction archive** or source checkout.
Install the pinned numerical dependencies with
`python -m pip install -r requirements.txt`. For graphical dependencies,
see the R setup below.

```sh
python scripts/reproduce.py --output reproduced
python scripts/reproduce.py --figures-only --output figures-from-archive
python -m unittest discover -s tests
```

The first command recomputes the benchmarks, diagnostic tables and all four R
figures with their 14 vector panels. The second renders the archived inputs.
Use `--skip-figures` for a numerical-only rerun. Each command requires a fresh
output directory, keeping the archived CSV and JSON files intact.

| Figure | Inputs |
| --- | --- |
| 1 | Prescribed geometry and contact-law parameters in the R source |
| 2 | Results from scripts 42, 45 and 48 |
| 3 | Results from script 51 |
| 4 | Results from scripts 45 and 51 and the fixed-order and isolated-breakpoint diagnostics |

The interface calls the same numerical rules as the original benchmarks.
The original standalone scripts and reference implementations remain available
for comparisons against the paper. Cite the software using `CITATION.cff`.

## Archived benchmark details

The CSV and JSON files are the saved results reported in the manuscript. Files 42, 45, 48 and 51 generate the inclined-wall, circular quadrature, circular closed-form and noncircular comparisons, respectively; 54 checks the noncircular boundary roots and circular limit. Protocols 41, 44, 47, 50 and 53 describe the corresponding runs. These synthetic configurations are not experimental observations.

Original results used NumPy 2.3.5, SciPy 1.16.3 and mpmath 1.3.0; the figure script used R 4.5, ggplot2 4.0.2, cowplot 1.2.0, svglite 2.2.2 and ragg 1.5.0. Output timings depend on machine and workload. A new run will overwrite saved CSV/JSON files: copy this directory before running the scripts. For the circular and noncircular angular benchmarks, high-precision references are rounded to binary64 before calculating the reported differences. The inclined-wall benchmark uses a separate segmented SciPy reference, with a 40-digit check at one prespecified height; the reference estimate is not a rigorous error bound.

Suggested order in the copied directory: python3 42_indicador_geometria.py; python3 45_indicador_gk.py; python3 48_forma_fechada.py; python3 51_teste_nao_circular.py; python3 54_verificar_limite_raizes.py. The Python scripts have fixed cases and write to their own directory. Run Rscript figures_r/render_contact_figures.R from this directory to regenerate the figures from the saved CSV files; only the input path has been adapted for this package.

The `review_checks` directory contains post hoc checks. To reproduce the inputs of Figure 4 after the primary scripts, run `python3 review_checks/isolated_breakpoints.py` and `python3 review_checks/fixed_segment_gauss.py` before the R script. The first writes both its JSON and CSV outputs; the second writes a CSV and a summary. The additional scripts `nonpenetrating_subset.py`, `check_decisive_references.py`, and `check_noncircular_maxima.py` reproduce the descriptive subset and precision checks. These analyses were added after inspection of the primary results and are exploratory.

Set up the Python dependencies with `python3 -m pip install -r requirements.txt` in a virtual environment. Install the R packages with `install.packages(c("ggplot2", "cowplot", "svglite", "ragg"))`. The plotting code uses Arial; if it is unavailable, install a licensed or system-provided Arial font before reproducing the publication typography. Cairo support is needed for the vector PDF device. Computed results can be checked without R; the listed R packages and fonts are needed only for figures. Timings are descriptive, are remeasured on each run, and will vary across machines and workloads.

Author: Renato Socodato, University of Porto. The code and data accompany the manuscript. The Python and R code in this repository is distributed under the MIT License (see `LICENSE`). The numerical data, figures and other non-code materials are not covered by that code licence; contact the author regarding their reuse.

The computed quantities describe idealized contact models and synthetic configurations; they are not measurements from biological specimens. The 52 circular and 18 noncircular cases, near-endpoint selection and post hoc checks are described in the manuscript and protocols. The strongest observed accuracy violations occur in strongly interpenetrating penalty-model configurations.

The 14 individual figure panels and four composite figures are in `figures_r/`; the numerical scripts save their output beside the scripts, so use a copy for re-execution. The scripts in `review_checks/` use the saved primary tables and include exploratory checks.
