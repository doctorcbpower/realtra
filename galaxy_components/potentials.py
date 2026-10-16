"""Analytic mass components, replacing the #ifdef blocks of the Fortran code.

Each component provides
    phi(R, z)  potential [(km/s)^2]            (arrays broadcast)
    vc2(R)     midplane circular speed squared [(km/s)^2]
so a model is simply a list of components whose contributions are summed.
"""
from dataclasses import dataclass

import numpy as np
from scipy import special

from .constants import GGDT, PI


class Component:
    def phi(self, R, z):
        raise NotImplementedError

    def vc2(self, R):
        raise NotImplementedError


@dataclass
class Hernquist(Component):
    """Hernquist sphere.  ``mass`` is the asymptotic total mass."""
    mass: float
    a: float

    def phi(self, R, z):
        return -GGDT * self.mass / (self.a + np.hypot(R, z))

    def menc(self, r):
        x = np.asarray(r) / self.a
        return self.mass * x * x / (1 + x) ** 2

    def vc2(self, R):
        return GGDT * self.menc(R) / R


@dataclass
class NFW(Component):
    """NFW halo with ``mass`` the mass inside r200 = c * rs."""
    mass: float
    rs: float
    c: float

    @property
    def _fofc(self):
        return np.log(1 + self.c) - self.c / (1 + self.c)

    def phi(self, R, z):
        x = np.hypot(R, z) / self.rs
        return -GGDT * self.mass / (self.rs * self._fofc) * np.log1p(x) / x

    def menc(self, r):
        x = np.asarray(r) / self.rs
        return self.mass / self._fofc * (np.log1p(x) - x / (1 + x))

    def vc2(self, R):
        return GGDT * self.menc(R) / R


@dataclass
class LogarithmicHalo(Component):
    v0: float
    rcore: float
    q: float = 1.0

    def phi(self, R, z):
        return 0.5 * self.v0 ** 2 * np.log(self.rcore ** 2 + R ** 2 + (z / self.q) ** 2)

    def vc2(self, R):
        return self.v0 ** 2 * R ** 2 / (self.rcore ** 2 + R ** 2)


def _expdisc_nodes(R, Rd, nglq=16):
    """Gauss-Legendre nodes/weights in k for the Hankel integral at radius R.

    Panels resolve (i) the log-spaced scale 1/Rd .. 1/z and (ii) the J0(kR)
    oscillations (period 2 pi / R), so the same rule is valid for all z.
    """
    kmax = 150.0 / Rd
    geom = np.geomspace(1e-4 / Rd, kmax, 60)
    nlin = int(np.ceil(kmax * R / PI)) + 1
    lin = np.linspace(0.0, kmax, nlin + 1) if nlin > 1 else np.array([0.0, kmax])
    edges = np.unique(np.concatenate(([0.0], geom, lin)))
    xg, wg = np.polynomial.legendre.leggauss(nglq)
    lo, hi = edges[:-1, None], edges[1:, None]
    k = (0.5 * (hi - lo) * xg + 0.5 * (hi + lo)).ravel()
    w = (0.5 * (hi - lo) * wg).ravel()
    return k, w


@dataclass
class ExponentialDisc(Component):
    """Razor-thin exponential disc, Sigma = M/(2 pi Rd^2) exp(-R/Rd).

    phi(R,z) = -G M int_0^inf J0(kR) exp(-k|z|) / (1+(k Rd)^2)^{3/2} dk
    """
    mass: float
    rd: float

    @property
    def sigma0(self):
        return self.mass / (2 * PI * self.rd ** 2)

    def phi(self, R, z):
        R = np.atleast_1d(np.asarray(R, dtype=float))
        z = np.abs(np.atleast_1d(np.asarray(z, dtype=float)))
        out = np.empty((R.size, z.size))
        for i, Ri in enumerate(R):
            k, w = _expdisc_nodes(Ri, self.rd)
            amp = w * special.j0(k * Ri) / (1 + (k * self.rd) ** 2) ** 1.5
            out[i] = -GGDT * self.mass * (np.exp(-np.outer(z, k)) @ amp)
        return out

    def vc2(self, R):
        """Freeman (1970): vc^2 = 4 pi G Sigma0 Rd y^2 [I0K0 - I1K1], y = R/2Rd."""
        y = 0.5 * np.asarray(R) / self.rd
        return (4 * PI * GGDT * self.sigma0 * self.rd * y * y *
                (special.i0(y) * special.k0(y) - special.i1(y) * special.k1(y)))
