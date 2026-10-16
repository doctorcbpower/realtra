from pathlib import Path

import numpy as np
from scipy import special

from galaxy_components import gas_disc, nbody
from galaxy_components.constants import GGDT, PI
from galaxy_components.params import read_params
from galaxy_components.potentials import ExponentialDisc

DATA = Path(__file__).parent / "data"


def test_exponential_disc_midplane_potential_matches_freeman():
    d = ExponentialDisc(mass=1.0, rd=2.0)
    R = np.array([0.1, 1.0, 4.0, 20.0, 150.0])
    phi = d.phi(R, np.array([1e-12]))[:, 0]
    y = R / (2 * d.rd)
    exact = -PI * GGDT * d.sigma0 * R * (special.i0(y) * special.k1(y) - special.i1(y) * special.k0(y))
    assert np.allclose(phi, exact, rtol=2e-4)


def test_exponential_disc_vc2_matches_potential_gradient():
    d = ExponentialDisc(mass=1.0, rd=2.0)
    R = np.array([1.0, 3.0, 8.0])
    h = 1e-3
    dphi = (d.phi(R + h, np.array([1e-12])) - d.phi(R - h, np.array([1e-12])))[:, 0] / (2 * h)
    assert np.allclose(R * dphi, d.vc2(R), rtol=1e-3)


def test_isolated_isothermal_sheet_is_sech2():
    # no background: rho = rho0 sech^2(z/z0), z0 = cs^2 / (pi G Sigma) = Sigma/(2 rho0)
    cs, sigma = 9.0, 0.01
    z = 10 ** np.linspace(-6, 1.5, 400)
    rho0, f, cum, _ = gas_disc.solve_ring(sigma, z, np.zeros_like(z), cs)
    z0 = cs ** 2 / (PI * GGDT * sigma)
    assert abs(sigma / (2 * rho0) / z0 - 1) < 2e-3
    assert np.allclose(f, 1 / np.cosh(z / z0) ** 2, rtol=5e-3, atol=1e-6)


def test_nbody_distribution_function_matches_fortran():
    p = read_params(DATA / "nbody_ref_params.txt")
    nb = nbody.build_nbody(p, verbose=False)
    ref = np.loadtxt(DATA / "nbody_ref_distribution_function.txt")
    n = len(ref)
    assert np.allclose(nb.etab[:n], ref[:, 0], rtol=1e-4)
    py = nb.ftab[0][:n] * nb.rscale_comp[0] ** 3
    m = ref[:, 1] > 0
    assert ((py > 0) == m).all()
    assert np.abs(py[m] / ref[m, 1] - 1).max() < 0.015  # ~0.9% observed


def test_vectorised_ring_solver_matches_scipy_reference():
    cs = 9.0
    z = 10 ** np.linspace(-15, 1.3, 500)
    # background: a mix of mild and strong vertical potentials
    R = np.array([0.5, 2.0, 8.0])
    phibz = np.stack([0.5 * k ** 2 * z ** 2 / (1 + z) for k in (30.0, 80.0, 150.0)])
    sig = np.array([0.02, 0.005, 0.0008])
    rho0, f, cum, _ = gas_disc.solve_rings(sig, z, phibz, cs)
    for i in range(3):
        r_ref, f_ref, cum_ref, _ = gas_disc.solve_ring(sig[i], z, phibz[i], cs)
        assert abs(rho0[i] / r_ref - 1) < 2e-4
        assert np.allclose(cum[i, -1], cum_ref[-1], rtol=2e-4)


def test_vectorised_solver_sech2():
    cs, sigma = 9.0, 0.01
    z = 10 ** np.linspace(-6, 1.5, 400)
    rho0, f, cum, _ = gas_disc.solve_rings(np.array([sigma]), z, np.zeros((1, z.size)), cs)
    z0 = cs ** 2 / (PI * GGDT * sigma)
    assert abs(sigma / (2 * rho0[0]) / z0 - 1) < 2e-3
