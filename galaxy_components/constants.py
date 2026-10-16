"""Physical constants and unit system.

Units follow the Fortran code: lengths in kpc, velocities in km/s and masses
in 1e10 Msol, so that G is the value below.  (The Fortran code also truncates
pi to 3.14159; here we use the exact value, a ~1e-6 relative difference.)
"""
import numpy as np

PI = np.pi

MPC_IN_M = 3.086e22
MSOL_IN_KG = 1.989e30
KMS_IN_MS = 1000.0
LUNIT = 1e-3  # kpc in Mpc
VUNIT = 1.0
MUNIT = 1e10
G_MKS = 6.67e-11

# G in kpc (km/s)^2 / (1e10 Msol)
GGDT = G_MKS * (MUNIT * MSOL_IN_KG) / (LUNIT * MPC_IN_M) / (VUNIT * KMS_IN_MS) ** 2

KB = 1.3806e-23  # Boltzmann constant, MKS
MP = 1.6726e-27  # proton mass, kg
GAMMA = 5.0 / 3.0

# Hydrogen / helium / metal mass fractions of the gas
XH, YHE, ZM = 1.0, 0.0, 0.0

RHOCRIT0 = 27.755e-9  # critical density, 1e10 Msol / kpc^3 / h^2
DELTAVIR = 200.0

LN10 = np.log(10.0)


def virial_radius(mvir, delta=DELTAVIR, rhoc=RHOCRIT0):
    return (3.0 * mvir / (4.0 * PI) / delta / rhoc) ** (1.0 / 3.0)


def virial_mass(vvir, delta=DELTAVIR, rhoc=RHOCRIT0):
    return np.sqrt(3.0 / (4.0 * PI) / delta / rhoc) * (vvir / np.sqrt(GGDT)) ** 3


def mean_molecular_weight(temperature, x=XH, y=YHE):
    """mu for neutral gas below 1e4 K, fully ionised above (as the Fortran code)."""
    mu_ion = 1.0 / (x / 1.0 + y / 4.0)
    if np.log10(temperature) <= 4:
        return mu_ion
    mu_e = 1.0 / (x / 1.0 + 2.0 * y / 4.0)
    return 1.0 / (1.0 / mu_ion + 1.0 / mu_e)


def sound_speed(temperature):
    """Isothermal sound speed in km/s."""
    mu = mean_molecular_weight(temperature)
    return np.sqrt(KB * temperature / mu / MP) * 1e-3
