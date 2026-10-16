"""Spherical N-body halo by Eddington inversion (isotropic, no disc in the potential).

    PYTHONPATH=. python examples/03_nbody_halo.py
"""
import numpy as np

from galaxy_components import nbody
from galaxy_components.presets import example_params

p = example_params(npart_comp=[30000, 0, 0, 0, 0])
nb = nbody.build_nbody(p, verbose=False)
comps = nbody.sample_nbody(p, nb, np.random.default_rng(3), verbose=False)
s = nbody.summarise(nb, comps)

pos = comps[0]["pos"].astype(float)
r = np.linalg.norm(pos, axis=1)
print(f"{len(pos)} particles of mass {comps[0]['mpart']:.3g}; {s['nvir']} inside r200 = {nb.r200:.1f} kpc")
print(f"velocity dispersion inside r200: {s['sigma']} km/s (total {s['sigma_tot']:.1f}); "
      f"radial dispersion {s['sigma_r']:.1f}, mean v_r {s['vr']:.2f}")
print("\n   r [kpc]    M_particles(<r)   M_model(<r)     (1e10 Msol)")
for rr in (5, 20, 80, 160):
    print(f"  {rr:7.0f}    {comps[0]['mpart'] * (r < rr).sum():12.3f}   {np.interp(rr, nb.rtab, nb.menc_tab[0]):12.3f}")
