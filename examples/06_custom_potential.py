"""Use the pieces on their own: compose potentials, and reuse the vertical-equilibrium solver
for a potential of your own (here a flattened logarithmic halo plus an exponential disc).

    PYTHONPATH=. python examples/06_custom_potential.py
"""
import numpy as np

from galaxy_components import gas_disc
from galaxy_components.constants import GGDT, PI, sound_speed
from galaxy_components.potentials import ExponentialDisc, Hernquist, LogarithmicHalo

# 1. compose analytic components and read off the circular-speed curve
parts = [Hernquist(mass=8.0, a=0.7), ExponentialDisc(mass=5.0, rd=3.0),
         LogarithmicHalo(v0=180.0, rcore=4.0, q=0.9)]
R = np.array([1.0, 3.0, 6.0, 12.0, 24.0])
print("R [kpc]       " + "  ".join(f"{R_:6.1f}" for R_ in R))
for c in parts:
    print(f"{type(c).__name__:<18}" + "  ".join(f"{np.sqrt(c.vc2(R_)):6.1f}" for R_ in R))
print(f"{'total':<18}" + "  ".join(f"{np.sqrt(sum(c.vc2(R_) for c in parts)):6.1f}" for R_ in R))

# 2. isothermal gas disc in that potential: any Phi(R, z) works, via solve_rings
z = 10 ** np.linspace(-6, 1.3, 300)                 # z grid, uniform in ln z
Rr = np.array([2.0, 6.0, 12.0])
phi = np.stack([sum(c.phi(np.array([r]), z)[0] if isinstance(c, ExponentialDisc) else c.phi(r, z)
                    for c in parts) for r in Rr])
phibz = phi - phi[:, :1]                            # potential relative to the midplane
cs = sound_speed(1e4)
sigma = 0.01 * np.exp(-Rr / 3.0)                    # imposed surface density [1e10 Msol/kpc^2]
rho0, f, cum, nev = gas_disc.solve_rings(sigma, z, phibz, cs)
print(f"\ngas disc (T = 1e4 K, solved in {nev} iterations)")
for r, s, r0, fr in zip(Rr, sigma, rho0, f):
    zhalf = z[np.argmax(fr < 0.5)]                  # height where density falls to half the midplane
    print(f"  R = {r:5.1f} kpc: midplane rho = {r0:.3e}, Sigma/2rho0 = {1e3 * s / 2 / r0:7.1f} pc, "
          f"half-density height = {1e3 * zhalf:7.1f} pc")
