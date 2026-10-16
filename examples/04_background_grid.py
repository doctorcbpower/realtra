"""AREPO background grid for ANY set of gas particles (no other module needed).

With --input, read gas from an existing GADGET HDF5 file; otherwise use a random gas ball.

    PYTHONPATH=. python examples/04_background_grid.py [--input snap.hdf5] [--out arepo_ic.hdf5]
"""
import argparse

import numpy as np

from galaxy_components import backgroundgrid
from galaxy_components.io import write_hdf5
from galaxy_components.system import read_system

ap = argparse.ArgumentParser()
ap.add_argument("--input")
ap.add_argument("--out")
args = ap.parse_args()

if args.input:
    gas = read_system(args.input).species[0]
    pos, vel, mass, u = gas["pos"], gas["vel"], gas["mass"], gas["ug"]
else:
    rng = np.random.default_rng(4)
    pos = rng.normal(size=(20000, 3)) * 30.0
    vel, mass, u = rng.normal(size=(20000, 3)) * 5.0, 1e-4, np.full(20000, 100.0)

out = backgroundgrid.add_background_grid(pos, vel, mass, u, boxsize_factor=5.0, nngb=64, verbose=False)
m = out["mass"]
print(f"{len(pos)} particles -> {len(m)} cells in a box of {out['boxsize']:.1f}")
print(f"cell side: {out['volume'].min() ** (1 / 3):.3g} .. {out['volume'].max() ** (1 / 3):.3g}")
print(f"mass in {np.sum(np.broadcast_to(mass, len(pos))):.6g} -> out {m.sum():.6g}; "
      f"thermal energy in {np.sum(np.broadcast_to(mass, len(pos)) * u):.6g} -> out {np.sum(m * out['u']):.6g}")
print("AREPO parameters:", out["suggest"])
if args.out:
    write_hdf5(args.out, {0: dict(pos=out["pos"], vel=out["vel"], mass=m, ug=out["u"])},
               boxsize=out["boxsize"], double=True)
