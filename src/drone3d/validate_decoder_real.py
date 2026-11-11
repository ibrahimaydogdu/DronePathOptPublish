"""Equivalence check (Review_Round_03 condition 3): on REAL legs, the minimum of the fast decoder over all
event permutations must equal the exact solver's optimum (same model, t_service = 0).

Usage: python -I validate_decoder_real.py <n_req> <instance_seed> [max_events]
Writes nothing; prints one line per instance. Results are recorded by hand in 04_Results/VALIDATION_LOG.md.
"""
import itertools
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
import numpy as np

from around import CityRaster
from energy import Hardware
from exact import solve_exact
from fast import FastProblem, build_table_sym
from geometry import load_city
from instances import DEPOT_A, load_pluto, make_instance
from legs import LegCost


def main():
    n_req, seed = int(sys.argv[1]), int(sys.argv[2])
    max_events = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    city = load_city(os.path.join(ROOT, "Evidence", "data", "NYC_Building_Footprints_raw.json"))
    pluto = load_pluto(os.path.join(ROOT, "Evidence", "data", "NYC_MapPLUTO_raw.json"))
    ras = CityRaster(city, res=5.0)
    hw = Hardware()
    inst = make_instance(city, pluto, DEPOT_A, n_req, seed)
    leg = LegCost(city, inst.nodes, model="F", raster=ras)
    E = build_table_sym(leg, len(inst.nodes), hw.q_max)
    fast = FastProblem(E, inst.requests, hw.q_max, hw.e_use_wh, kmax=24)
    if fast.D > max_events:
        print(f"seed {seed}: {fast.D} events > {max_events}, skipped")
        return
    opt, trips = solve_exact(leg, inst.requests, hw.q_max, hw.e_use_wh)
    best = float("inf")
    t = time.time()
    for perm in itertools.permutations(range(fast.D)):
        x = np.empty(fast.D)
        for rank, k in enumerate(perm):
            x[k] = (rank + 0.5) / fast.D
        e, ok, _ = fast.feasible_energy(x)
        if ok and e < best:
            best = e
    print(f"seed {seed} n_req {n_req} events {fast.D} exact {opt:.6f} decoder-enumeration {best:.6f} "
          f"diff {abs(opt - best):.2e} trips {len(trips)} ({time.time() - t:.0f}s)")


if __name__ == "__main__":
    main()
