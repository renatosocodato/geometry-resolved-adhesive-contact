"""Recompute the paper's tables and figures in a new directory."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = (
    "42_indicador_geometria.py", "45_indicador_gk.py", "48_forma_fechada.py",
    "51_teste_nao_circular.py", "54_verificar_limite_raizes.py",
)
CHECKS = (
    "isolated_breakpoints.py", "fixed_segment_gauss.py", "nonpenetrating_subset.py",
    "check_decisive_references.py", "check_noncircular_maxima.py",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new directory")
    parser.add_argument("--figures-only", action="store_true", help="use the archived CSV inputs")
    parser.add_argument("--skip-figures", action="store_true", help="run numerical scripts only")
    args = parser.parse_args()
    if args.figures_only and args.skip_figures:
        parser.error("choose at most one figure option")
    target = args.output.resolve()
    if target.exists():
        parser.error("output directory already exists")
    if not args.skip_figures and shutil.which("Rscript") is None:
        parser.error("Rscript is required for figure rendering")

    target.mkdir(parents=True)
    for relative in ("review_checks", "figures_r"):
        (target / relative).mkdir()
    if args.figures_only:
        for pattern in ("*.csv", "review_checks/*.csv"):
            for source in ROOT.glob(pattern):
                shutil.copy2(source, target / source.relative_to(ROOT))
    else:
        for relative in (*BENCHMARKS, *(f"review_checks/{name}" for name in CHECKS)):
            shutil.copy2(ROOT / relative, target / relative)
            subprocess.run([sys.executable, relative], cwd=target, check=True)
    if not args.skip_figures:
        relative = "figures_r/render_contact_figures.R"
        shutil.copy2(ROOT / relative, target / relative)
        subprocess.run(["Rscript", "--vanilla", relative], cwd=target, check=True)


if __name__ == "__main__":
    main()
