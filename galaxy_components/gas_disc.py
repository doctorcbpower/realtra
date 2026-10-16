"""Isothermal gas disc in hydrostatic equilibrium (port of mk_gas_disc.f90).

At each radius R the vertical structure solves

    rho(z) = rho0 exp(-(Phi_bg(R,z) - Phi_bg(R,0) + Phi_g(z)) / cs^2)
    Phi_g'' = 4 pi G rho(z),   Phi_g(0) = Phi_g'(0) = 0

(the gas's own vertical self-gravity, in the local slab approximation), with
rho0 chosen so that 2 * int_0^inf rho dz = Sigma(R).  Sigma(R) is imposed.
"""
from dataclasses import dataclass

import numpy as np
from scipy import integrate, optimize, special

from .constants import GAMMA, GGDT, PI, sound_speed
from .model import GalaxyModel
from .params import Params
from .potentials import ExponentialDisc


def sigma_shape(R, rhole, rref, rd):
    """Sigma(R)/sigma0 : exponential, optionally with an inner hole."""
    R = np.asarray(R, dtype=float)
    hole = R * R / (R * R + rhole * rhole) if rhole > 0 else 1.0
    return hole * np.exp(-(R - rref) / rd)


@dataclass
class GasDisc:
    model: GalaxyModel
    cs: float
    sigma0: float
    rtab: np.ndarray
    ztab: np.ndarray
    pot: np.ndarray       # background potential, (nr, nz)
    rho0: np.ndarray      # midplane density (nr,)
    lrho: np.ndarray      # log10 rho(R,z), (nr, nz)
    mztab: np.ndarray     # normalised cumulative vertical mass, (nr, nz)
    mrtab: np.ndarray     # cumulative gas mass inside R (nr,)
    vrot: np.ndarray      # rotation speed (nr,)
    vcomp: dict           # components of vrot^2
    niter: np.ndarray     # function evaluations per ring

    @property
    def gas_mass(self):
        return self.mrtab[-1]


def _vertical_solve(rho0, zgrid, phibz, cs):
    """Integrate the vertical Poisson equation on zgrid; returns Phi_g(z)."""
    t0, t1 = np.log(zgrid[0]), np.log(zgrid[-1])
    c4 = 4 * PI * GGDT * rho0

    def rhs(t, y):
        z = np.exp(t)
        return [z * y[1], z * c4 * np.exp(-(np.interp(z, zgrid, phibz) + y[0]) / cs ** 2)]

    sol = integrate.solve_ivp(rhs, (t0, t1), [0.0, 0.0], method="RK45",
                              rtol=1e-9, atol=1e-12, t_eval=np.log(zgrid))
    if not sol.success or sol.y.shape[1] != zgrid.size:
        raise RuntimeError(f"vertical integration failed: {sol.message}")
    return sol.y[0]


def _vertical_integral(zgrid, phibz, phig, cs):
    """f(z) = exp(-(Phi_bg+Phi_g)/cs^2) and its cumulative integral over z."""
    f = np.exp(-(phibz + phig) / cs ** 2)
    # integrate in ln z (f*z is smooth there); [0, z0] contributes ~ f0*z0
    cum = integrate.cumulative_simpson(f * zgrid, x=np.log(zgrid), initial=0.0) + f[0] * zgrid[0]
    return f, cum


def solve_ring(sigma_R, zgrid, phibz, cs, tol=1e-8):
    """Find rho0 s.t. 2*int rho dz = sigma_R.  Returns (rho0, f, cum, n_evals)."""
    nev = 0

    def total(rho0):
        nonlocal nev
        nev += 1
        phig = _vertical_solve(rho0, zgrid, phibz, cs)
        f, cum = _vertical_integral(zgrid, phibz, phig, cs)
        return f, cum

    f0, cum0 = _vertical_integral(zgrid, phibz, np.zeros_like(zgrid), cs)
    guess = sigma_R / (2 * cum0[-1])

    def resid(lrho):
        f, cum = total(np.exp(lrho))
        return np.log(2 * np.exp(lrho) * cum[-1] / sigma_R)

    l0 = np.log(guess)
    # self-gravity only compresses the disc, so rho0 >= guess in practice
    a, b = l0 - 0.1, l0 + 2.0
    fa, fb = resid(a), resid(b)
    while fa * fb > 0:
        a, b = a - 1.0, b + 1.0
        fa, fb = resid(a), resid(b)
        if b - l0 > 30:
            raise RuntimeError("could not bracket rho0")
    lrho = optimize.brentq(resid, a, b, xtol=tol, rtol=tol)
    rho0 = np.exp(lrho)
    f, cum = total(rho0)
    return rho0, f, cum, nev


def build_gas_disc(p: Params, model: GalaxyModel, verbose=True) -> GasDisc:
    log = print if verbose else (lambda *a, **k: None)
    rd = model.rdisc
    cs = sound_speed(p.temp)
    log(f"Sound speed [km/s]: {cs:.4f}")

    rtab = 10 ** np.linspace(np.log10(rd) - 3.0, np.log10(rd) + 2.0, p.nrbin)
    ztab = 10 ** np.linspace(p.log_zmin, p.log_zmax, p.nzbin)

    # Surface density normalisation: gas mass = mdisc - mstar (or sigma0 if given)
    shape = lambda R: sigma_shape(R, model.rdisc_hole, model.rdisc_ref, rd)
    rmin = max(rtab[0], model.rdisc_hole)
    lr = np.linspace(np.log(rmin), np.log(rtab[-1]), 20001)
    integral = integrate.trapezoid(2 * PI * np.exp(2 * lr) * shape(np.exp(lr)), lr)
    sigma0 = p.sigma0 if p.sigma0 not in (0.0, 1.0) else (model.mdisc - model.mstar) / integral
    log(f"Sigma0 (gas) [Msol/pc^2]: {sigma0 * 1e4:.4f}; gas mass [1e10 Msol]: {sigma0 * integral:.6g}")

    log("Computing background potential...")
    pot = model.phi_background(rtab, ztab)

    nr, nz = pot.shape
    rho0 = np.empty(nr)
    lrho = np.empty((nr, nz))
    mztab = np.empty((nr, nz))
    niter = np.empty(nr, dtype=int)
    for i in range(nr):
        phibz = pot[i] - pot[i, 0]
        sig = sigma0 * shape(rtab[i])
        r0, f, cum, nev = solve_ring(sig, ztab, phibz, cs)
        rho0[i], niter[i] = r0, nev
        lrho[i] = np.log10(np.maximum(r0 * f, 1e-300))
        m = cum / cum[-1]
        m[1.0 - m < 5e-3] = 1.0
        mztab[i] = m
        log(f"ring {i + 1:4d} R={rtab[i]:.5e} rho0={r0:.5e} Sigma/2rho0={sig / 2 / r0:.4e} evals={nev}")

    # cumulative gas mass vs R (honours the hole, unlike the Fortran analytic formula)
    mrtab = integrate.cumulative_trapezoid(2 * PI * sigma0 * rtab * shape(rtab), rtab, initial=0.0)

    vcomp = rotation_curve(p, model, rtab, lrho, sigma0, cs)
    vrot = np.sqrt(np.maximum(vcomp["total"], 0.0))
    return GasDisc(model, cs, sigma0, rtab, ztab, pot, rho0, lrho, mztab, mrtab, vrot, vcomp, niter)


def rotation_curve(p, model, rtab, lrho, sigma0, cs):
    """vrot^2 = vc^2(background) + vc^2(gas) + cs^2 dln(rho_mid)/dln(R)."""
    gas = ExponentialDisc(sigma0 * 2 * PI * model.rdisc ** 2, model.rdisc)
    if model.rdisc_hole > 0:
        print("Warning: gas-disc rotation uses the Freeman exponential formula, ignoring the hole")
    v_bg = model.vc2_background(rtab)
    v_gas = gas.vc2(rtab)
    lr = np.log10(rtab)
    vp = cs ** 2 * np.gradient(lrho[:, 0], lr)  # d log10 rho / d log10 R  (base cancels)
    total = v_bg + v_gas + vp
    return {"background": v_bg, "gas": v_gas, "pressure": vp, "total": total}


def sample_particles(p: Params, disc: GasDisc, ndisc, rng):
    """Draw gas particles; returns pos, vel (float32, N x 3), u, mpart."""
    model, rtab, ztab = disc.model, disc.rtab, disc.ztab
    rd = model.rdisc

    # outer edge of the disc is 25 scale lengths (as the Fortran code)
    rcut = min(rtab[-1], 10 ** (np.log10(rd) + 1.4))
    mkeep = np.interp(rcut, rtab, disc.mrtab)
    r0 = p.fhole * rd
    m0 = np.interp(r0, rtab, disc.mrtab) if r0 > 0 else 0.0
    mpart = (mkeep - m0) / ndisc
    print(f"Disc mass to R={rcut:.3f} kpc: {mkeep:.6g}; particle mass {mpart:.4e}")

    if not p.ispoisson:
        raise NotImplementedError("glass-file sampling (ispoisson=.false.) is not ported yet")

    # radii by inversion of M(<R), rejecting anything inside the hole
    r = np.empty(0)
    row = np.empty(0, dtype=int)
    while r.size < ndisc:
        need = int((ndisc - r.size) * 1.2) + 10
        mt = rng.random(need) * mkeep
        rr = np.interp(mt, disc.mrtab, rtab)
        rr = rr[rr >= r0]
        r = np.concatenate((r, rr))
    r = r[:ndisc]
    row = np.clip(np.searchsorted(rtab, r, side="right") - 1, 0, rtab.size - 1)

    # heights by inversion of the normalised vertical mass profile of the ring below
    u = rng.random(ndisc)
    zabs = np.empty(ndisc)
    for i in np.unique(row):
        sel = row == i
        m = disc.mztab[i]
        if p.snap_z:
            j = np.minimum(np.searchsorted(m, u[sel], side="right"), ztab.size - 1)
            zabs[sel] = ztab[j]
        else:
            zabs[sel] = np.exp(np.interp(u[sel], m, np.log(ztab)))
    z = np.where(rng.random(ndisc) < 0.5, -zabs, zabs)
    phi = 2 * PI * rng.random(ndisc)

    vc = np.interp(r, rtab, disc.vrot)
    pos = np.column_stack((r * np.cos(phi), r * np.sin(phi), z))
    vel = np.column_stack((-vc * np.sin(phi), vc * np.cos(phi), np.zeros(ndisc)))
    ug = np.full(ndisc, disc.cs ** 2 / (GAMMA - 1))
    return pos.astype(np.float32), vel.astype(np.float32), ug.astype(np.float32), mpart
