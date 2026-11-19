"""Run one calculation and emit JSON."""

import argparse
from dataclasses import asdict
import json

from . import contact_force


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height-um", type=float, required=True)
    parser.add_argument("--radius-um", type=float)
    parser.add_argument("--rx-um", type=float)
    parser.add_argument("--ry-um", type=float)
    parser.add_argument("--q-um", type=float, default=0.0)
    parser.add_argument("--length-um", type=float, default=0.75)
    parser.add_argument("--stiffness-pn-um3", dest="stiffness_pn_um3", type=float, default=20.0)
    parser.add_argument("--adhesion-pn-um", type=float, default=0.002)
    parser.add_argument("--range-um", type=float, default=0.1)
    parser.add_argument("--tolerance-pn", type=float, default=1e-8)
    parser.add_argument("--method", choices=("gk", "guarded", "all"), default="all")
    try:
        result = contact_force(**vars(parser.parse_args()))
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(asdict(result), indent=2, allow_nan=False))
