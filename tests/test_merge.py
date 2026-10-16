from pathlib import Path

import h5py
import numpy as np
import pytest

from galaxy_components import merge
from galaxy_components.constants import GGDT
from galaxy_components.system import System, read_system

DATA = Path(__file__).parent / "data"


def _toy(name, n=200, seed=0, gas=True):
    rng = np.random.default_rng(seed)
    sp = {1: dict(pos=rng.normal(size=(n, 3)) * 10, vel=rng.normal(size=(n, 3)) * 50, mass=0.5)}
    if gas:
        sp[0] = dict(pos=rng.normal(size=(n, 3)) * 5, vel=rng.normal(size=(n, 3)),
                     mass=rng.random(n) + 0.1, ug=rng.random(n) + 1)
    return System(name, sp, {"m200": 100.0})


def test_rotation_and_place():
    R = merge.rotation_matrix([30, 50, 70])
    assert np.allclose(R @ R.T, np.eye(3)) and np.isclose(np.linalg.det(R), 1)
    s = _toy("a")
    p = merge.place(s, position=[10, 0, 0], velocity=[0, 5, 0], euler=[20, 40, 60])
    r0 = np.linalg.norm(s.species[1]["pos"], axis=1)
    r1 = np.linalg.norm(p.species[1]["pos"] - [10, 0, 0], axis=1)
    assert np.allclose(r0, r1)
    assert np.allclose(p.species[1]["vel"].mean(0) - [0, 5, 0],
                       merge.rotation_matrix([20, 40, 60]) @ s.species[1]["vel"].mean(0))


@pytest.mark.parametrize("ecc", [1.0, 0.7, 1.5])
def test_two_body_orbit(ecc):
    m1, m2, r0, q = 100.0, 40.0, 120.0, 25.0   # r0 < apocentre (141) for e=0.7
    x1, v1, x2, v2 = merge.orbit_two_body(m1, m2, r0, q, ecc)
    assert np.allclose(m1 * x1 + m2 * x2, 0) and np.allclose(m1 * v1 + m2 * v2, 0)
    rel, vel = x2 - x1, v2 - v1
    mu = GGDT * (m1 + m2)
    p = q * (1 + ecc)
    assert np.isclose(np.linalg.norm(rel), r0)
    assert np.isclose(0.5 * vel @ vel - mu / r0, mu * (ecc ** 2 - 1) / (2 * p))   # energy of the conic
    assert np.isclose(np.linalg.norm(np.cross(rel, vel)), np.sqrt(mu * p))        # angular momentum
    assert rel @ vel < 0                                                           # incoming


def test_stitch_conserves_and_maps_types():
    a, b = _toy("a", seed=1), _toy("b", seed=2)
    out = merge.stitch([merge.place(a, [50, 0, 0]), merge.place(b, [-50, 0, 0])],
                       type_maps=[None, {1: 2}], verbose=False)
    assert set(out) == {0, 1, 2}
    assert len(out[0]["pos"]) == 400 and len(out[1]["pos"]) == 200 and len(out[2]["pos"]) == 200
    assert np.isclose(np.sum(out[0]["mass"]), a.mass(0) + b.mass(0))
    assert out[1]["mass"] == 0.5                         # uniform mass collapses to a scalar


def test_merge_comps_options():
    s = _toy("a", n=4000)
    rng = np.random.default_rng(0)
    t = merge.taper_profile(s, 0.8, rng)
    assert len(t.species[1]["pos"]) < len(s.species[1]["pos"]) and "taper_radius" in t.meta
    assert merge.excise_types(s, [1]).species.keys() == {0}
    assert np.isclose(merge.reduce_stellar_mass(s, 1, 0.5).species[1]["mass"], 0.25)


def test_roundtrip_and_background_grid(tmp_path):
    a, b = _toy("a", seed=3), _toy("b", seed=4)
    placed = [merge.place(a, [40, 0, 0]), merge.place(b, [-40, 0, 0])]
    species = merge.stitch(placed, verbose=False)
    gas_mass = np.sum(species[0]["mass"])
    species, box = merge.finalise(species, background_grid=True, verbose=False)
    path = tmp_path / "m.hdf5"
    merge.write_merged(path, species, box, verbose=False)
    with h5py.File(path) as f:
        assert f["Header"].attrs["BoxSize"] == pytest.approx(box) and box > 0
        assert np.isclose(f["PartType0/Masses"][:].sum(), gas_mass, rtol=1e-6)
        ids = np.concatenate([f[f"PartType{t}/ParticleIDs"][:] for t in (0, 1)])
        assert len(np.unique(ids)) == len(ids)
        for t in (0, 1):
            assert np.all((f[f"PartType{t}/Coordinates"][:] > 0) & (f[f"PartType{t}/Coordinates"][:] < box))
    back = read_system(path)
    assert np.isclose(back.mass(0), gas_mass, rtol=1e-6)


def test_merge_from_config_end_to_end(tmp_path):
    (tmp_path / "gal.txt").write_text((DATA / "galaxy_params.txt").read_text())
    (tmp_path / "m.toml").write_text(f"""
output = "merger.hdf5"
seed = 5
[background_grid]
enabled = true
[orbit]
systems = [0, 1]
r_init = 400.0
r_peri = 30.0
ecc = 1.0
[[system]]
name = "g1"
params = "gal.txt"
euler = [0, 30, 0]
[system.overrides]
ndisc = 400
nhalo = 600
npart_comp = [500, 0, 0, 0, 0]
[[system]]
name = "g2"
params = "gal.txt"
type_map = {{ 1 = 2 }}
[system.options]
halo_core = 0.1
[system.overrides]
ndisc = 300
nhalo = 400
npart_comp = [300, 0, 0, 0, 0]
""")
    out = merge.merge_from_config(tmp_path / "m.toml", tmp_path / "out", verbose=False)
    with h5py.File(out) as f:
        n = f["Header"].attrs["NumPart_Total"]
        assert n[1] > 0 and n[2] > 0 and n[0] > 0 and f["Header"].attrs["BoxSize"] > 0
        assert f["Header"].attrs["Flag_DoublePrecision"] == 1


def _sphere(name, centre, u, n=6000, radius=40.0, seed=0, mass_total=1.0):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 3))
    x = x / np.linalg.norm(x, axis=1)[:, None] * radius * rng.random(n)[:, None] ** (1 / 3)
    return System(name, {0: dict(pos=x + centre, vel=np.zeros((n, 3)), mass=mass_total / n,
                                 ug=np.full(n, float(u)))}, {})


def test_gas_overlap_keeps_higher_pressure_gas():
    hot = _sphere("hot", [0, 0, 0], u=4.0, seed=1)
    cold = _sphere("cold", [50, 0, 0], u=1.0, seed=2)          # spheres of radius 40 overlap in x in [10, 40]
    (h2, c2), stats = merge.resolve_gas_overlap([hot, cold], verbose=False)
    # the cold gas is cut where it lies inside the hot sphere's bulk ...
    assert stats["cold"]["removed"] > 0.15 * stats["cold"]["of"]
    xc = c2.species[0]["pos"]
    inside = np.linalg.norm(xc, axis=1) < 30                      # well inside the hot sphere
    assert inside.sum() < 0.02 * (np.linalg.norm(cold.species[0]["pos"], axis=1) < 30).sum() + 5
    # ... while the hot gas is essentially untouched
    assert stats["hot"]["removed"] < 0.02 * stats["hot"]["of"]
    # swapping the roles swaps the outcome
    hot2 = _sphere("hot", [50, 0, 0], u=4.0, seed=1)
    cold2 = _sphere("cold", [0, 0, 0], u=1.0, seed=2)
    _, st2 = merge.resolve_gas_overlap([cold2, hot2], verbose=False)
    assert st2["cold"]["removed"] > 0.15 * st2["cold"]["of"] and st2["hot"]["removed"] < 0.02 * st2["hot"]["of"]


def test_gas_overlap_leaves_separate_systems_alone():
    a = _sphere("a", [0, 0, 0], u=3.0, seed=3)
    b = _sphere("b", [500, 0, 0], u=1.0, seed=4)
    _, stats = merge.resolve_gas_overlap([a, b], verbose=False)
    assert stats["a"]["removed"] == 0 and stats["b"]["removed"] == 0
