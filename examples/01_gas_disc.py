"""Gas disc in hydrostatic equilibrium: potential -> vertical structure -> rotation curve -> particles.

    PYTHONPATH=. python examples/01_gas_disc.py [--out gas_disc.hdf5]
"""
import argparse

import numpy as np

from galaxy_components import gas_disc
from galaxy_components.io import write_hdf5
from galaxy_components.model import build_model
from galaxy_components.presets import example_params

ap = argparse.ArgumentParser()
ap.add_argument("--out", help="optionally write the particles to a GADGET HDF5 file")
args = ap.parse_args()

p = example_params(ndisc=4000)
model = build_model(p, verbose=False)          # halo/bulge/stellar-disc components + derived radii
disc = gas_disc.build_gas_disc(p, model, verbose=False)

print(f"R200 = {model.r200:.1f} kpc, disc scale length = {model.rdisc:.2f} kpc, "
      f"gas mass = {disc.gas_mass:.3f} (1e10 Msol), cs = {disc.cs:.2f} km/s")
print("\n   R [kpc]   vc_background  vc_gas   vrot   [km/s]     scale height [pc]")
for R in (0.5, 1, 2, 4, 8):
    i = np.argmin(abs(disc.rtab - R))
    sigma = disc.sigma0 * np.exp(-disc.rtab[i] / model.rdisc)
    h = sigma / (2 * disc.rho0[i])             # Sigma = 2 rho0 h   (sech^2-like sheet)
    print(f"  {disc.rtab[i]:7.2f}   {np.sqrt(disc.vcomp['background'][i]):9.1f}   "
          f"{np.sqrt(disc.vcomp['gas'][i]):6.1f}  {disc.vrot[i]:6.1f}      {1e3 * h:10.1f}")

pos, vel, ug, mpart = gas_disc.sample_particles(p, disc, p.ndisc, np.random.default_rng(1))
R = np.hypot(pos[:, 0], pos[:, 1])
print(f"\nsampled {len(pos)} particles of mass {mpart:.3g}; median R = {np.median(R):.2f} kpc, "
      f"rms z = {np.sqrt(np.mean(pos[:, 2] ** 2)):.3f} kpc")
if args.out:
    write_hdf5(args.out, {0: dict(pos=pos, vel=vel, ug=ug, mass=mpart)})
