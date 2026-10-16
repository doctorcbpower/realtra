"""Build two galaxies in memory, put them on an orbit, resolve the gas interface and stitch.

    PYTHONPATH=. python examples/05_stitch_systems.py [--out merger.hdf5]
"""
import argparse

from galaxy_components import merge
from galaxy_components.presets import example_params
from galaxy_components.system import build_system

ap = argparse.ArgumentParser()
ap.add_argument("--out", help="write the stitched IC (with an AREPO background grid)")
args = ap.parse_args()

small = dict(ndisc=1500, nhalo=3000, npart_comp=[6000, 0, 0, 0, 0])
a = build_system(example_params(**small), "primary", seed=1, halo_core=0.1)
b = build_system(example_params(m200=30.0, **small), "secondary", seed=2, halo_core=0.1)
print({s.name: (s.counts(), round(s.mass(), 2)) for s in (a, b)})

x1, v1, x2, v2 = merge.orbit_two_body(a.meta["m200"], b.meta["m200"], r_init=350.0, r_peri=30.0, ecc=1.0)
placed = [merge.place(a, x1, v1), merge.place(b, x2, v2, euler=(0, 60, 0))]   # tilt the second disc
placed, stats = merge.resolve_gas_overlap(placed, verbose=True)                 # pressure-balance the interface
species = merge.stitch(placed, type_maps=[None, {1: 2}], verbose=False)         # keep the 2nd halo as PartType2
species, boxsize = merge.finalise(species, background_grid=bool(args.out), verbose=False)
print({t: len(s["pos"]) for t, s in species.items()}, f"box {boxsize:.0f}" if boxsize else "")
if args.out:
    merge.write_merged(args.out, species, boxsize)
