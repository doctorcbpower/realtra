"""Hot gas halo: Komatsu & Seljak hydrostatic solution, entropy core and zero-pressure edge.

    PYTHONPATH=. python examples/02_gas_halo.py
"""
import numpy as np

from galaxy_components import gas_halo
from galaxy_components.presets import example_params

p = example_params()
ks = gas_halo.build_gas_halo(p, verbose=False)                      # same dark halo as the disc sees
cored = gas_halo.add_entropy_core(p, ks, 0.1, verbose=False)        # entropy floor -> density core
edged = gas_halo.add_entropy_core(p, ks, 0.1, verbose=False, zero_edge=True)   # P(rtrunc r200) = 0

mu_mp = 0.61 * 1.6726e-24
to_cgs = 1.989e43 / 3.086e21 ** 3                                   # 1e10 Msol/kpc^3 -> g/cm^3
print(f"K&S: gamma = {ks.gamma_poly:.3f}, central T = {ks.temp0:.3g} K, gas mass in r200 = {ks.mgas:.1f}")
print("\n   r [kpc]    n [cm^-3]  (K&S / cored / cored+edge)     T [K]  (K&S / cored / cored+edge)")
for r in (1.0, 10.0, 50.0, 150.0, 220.0):
    n = [h.rho(r) * to_cgs / mu_mp for h in (ks, cored, edged)]
    T = [h.u(r) * gas_halo.U_TO_TEMP * gas_halo.MU_HALO for h in (ks, cored, edged)]
    print(f"  {r:7.1f}    " + "  ".join(f"{v:9.2e}" for v in n) + "      " + "  ".join(f"{v:9.2e}" for v in T))

pos, vel, ug, mpart = gas_halo.sample_gas_halo(p, edged, np.random.default_rng(2), verbose=False)
print(f"\nsampled {len(pos)} particles (mass {mpart:.3g}); outermost at "
      f"{np.linalg.norm(pos, axis=1).max():.0f} kpc = {p.rtrunc:.1f} r200")
