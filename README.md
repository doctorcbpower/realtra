# build_galaxy_components

Equilibrium initial conditions for galaxy simulations: a **gas disc**, a hot **gas halo** and spherical
**N-body halo/bulge** components, built consistently in one dark-matter potential, stitched together
into one or more **systems** (e.g. for a merger), and optionally converted to an **AREPO** initial condition with a
background grid. Developed by [A/Prof Chris Power (ICRAR/UWA)](https://www.icrar.org/people/cpower/).

The maintained implementation is the Python package [`galaxy_components/`](galaxy_components/). The original
Fortran codes it grew out of are kept in [`legacy/`](legacy/) (with their own README).

```
params file ──► disc     gas disc in hydrostatic equilibrium (analytic potential)   ┐
 (namelists)   ──► halo     hot gas halo, hydrostatic in the same dark halo            ├─► system ─► merge ─► IC file
               ──► nbody    spherical isotropic halo/bulge, Eddington inversion        ┘   (compounded)  (orbit, overlap,
                                                                                                      background grid)
```

## Contents
[Install](#install) · [Quick start](#quick-start) · [Parameter file](#parameter-file) · [Commands](#commands) ·
[What each component models](#what-each-component-models) · [Merging systems](#merging-systems) ·
[Use the modules directly](#use-the-modules-directly) · [Tests](#tests) · [Differences from the Fortran](#differences-from-the-fortran) ·
[Limitations](#limitations) · [References](#references)

## Install

Python ≥ 3.11. From a clone of this repository:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[test]"        # numpy>=2, scipy>=1.12, h5py, f90nml (+ pytest); editable, so edits take effect
pytest -q                       # 49 tests, ~20 s
galaxy-components --help        # console script (also: python -m galaxy_components ...)
```
`pip install .` (non-editable) works too, and `pip wheel .` builds a wheel containing only the package. After installing, the scripts in
[`examples/`](examples/) run from anywhere (without installing, run them from the repository root with `PYTHONPATH=.`).

Units throughout: **kpc, km/s, 10¹⁰ M☉** (so G = 4.3×10⁴ kpc (km/s)² / 10¹⁰ M☉).

## Quick start

```bash
# a whole galaxy (gas disc + gas halo + N-body halo), with an AREPO background grid
galaxy-components galaxy examples/galaxy_params.txt --outdir out --background-grid --halo-core 0.1 --zero-edge

# a two-galaxy parabolic merger from a TOML description
galaxy-components merge examples/merger.toml --outdir out

# or script the pieces yourself
python examples/02_gas_halo.py
```
[`examples/galaxy_params.txt`](examples/galaxy_params.txt) is a ready-made parameter file for a 10¹² M☉ galaxy.

## Parameter file

A Fortran-style namelist file (see [`examples/galaxy_params.txt`](examples/galaxy_params.txt)) is read unchanged: the same file format
is used by the legacy Fortran programs. Groups and the fields used:

| group | fields | meaning |
|---|---|---|
| `&HALO` | `m200`, `v200`, `c200`, `qflat`, `rcore` | virial mass [10¹⁰ M☉] (or `v200`), concentration; log-halo flattening/core |
| `&DISC` | `mdisc`, `mstar`, `fdisc`, `fhole`, `fdref`, `sigma0` | disc mass (>0: fraction of M200; ≤0: absolute), stellar fraction, scale length as a fraction of the halo scale, inner hole, … |
| `&BULGE` | `mbulge`, `fbulge` | Hernquist bulge mass and scale |
| `&GAS_DISC` | `temp`, `ndisc` | gas-disc temperature [K], particle count |
| `&GAS_HALO` | `fbaryon`, `nhalo`, `rtrunc`, `form`, `v_radial` | baryon fraction, particles inside r200, sampling edge / r200, `'nfw'` or `'hernquist'` |
| `&BLACK_HOLE` | `mbh` | central black-hole mass |
| `&NBODY` | `ncomp`, `npart_comp`, `fmass_comp`, `form_comp`, `rscale_comp`, `rtrunc_comp`, `istaper` | N-body components (`'nfw'`, `'hernquist'`, `'burkert'`) |
| `&OUTPUT` | `disc_file`, `halo_file`, `nbody_file`, `snap_format` | output names; **only HDF5 (`snap_format=3`) is written** |
| `&MISC` | `ispoisson` | must be `.true.` (glass-file sampling is not ported) |
| `&MODEL` | `halo_type`, `galic`, `halo_core`, `halo_zero_edge`, `nrbin`, `nzbin`, `snap_z` | **Python only** (ignored by the Fortran): replaces the Fortran `-D` compile flags and holds new options |

`&MODEL`: `halo_type = 'hernquist' | 'nfw' | 'logarithmic'` (default `'hernquist'`), `galic = .true.` (GalIC
conventions: the Hernquist halo has total mass `mhalo` and scale `a = r200/c200 · √(2 f(c))`), `halo_core` (entropy
floor, below), `halo_zero_edge` (zero-pressure edge, below), `nrbin`/`nzbin` (gas-disc mesh, default 150 × 500).

## Commands

`galaxy-components <command> PARAMS [--outdir DIR] [--seed N]` (or `python -m galaxy_components ...`)

| command | builds | notable options |
|---|---|---|
| `disc` | gas disc (+ `vrot_components.txt`, `disc_parameters.txt`) | |
| `nbody` | spherical N-body components (+ `distribution_function.txt`) | `--fix-legacy` |
| `halo` | gas halo | `--halo-core F`, `--zero-edge`, `--legacy`, `--cold` |
| `galaxy` / `system` | disc + gas halo + N-body in one file | all `halo` options, plus `--naive-merge`, `--no-nbody`, `--background-grid`, `--boxsize-factor`, `--nngb`, `--name` |
| `merge CONFIG.toml` | several systems stitched into one file | see [Merging systems](#merging-systems) |

## What each component models

**Gas disc** (`disc`). Isothermal gas in hydrostatic equilibrium in the analytic potential of the halo, optional bulge
and stellar exponential disc. At each radius the vertical structure solves
ρ(z) = ρ₀ exp(−(Φ_bg(R,z) − Φ_bg(R,0) + Φ_g(z))/c_s²), Φ_g'' = 4πGρ, with ρ₀ set so the surface density follows the imposed
exponential. The gas's own gravity enters only through that local vertical term and, for the rotation curve, the Freeman
formula; the rotation curve is v_rot² = v_c²(background) + v_c²(gas) + c_s² dlnρ_mid/dlnR. All rings are solved together
with a vectorised RK4 (≈ 2 s for the default mesh). The stellar-disc potential is a Hankel-transform integral.

**N-body components** (`nbody`). Spherical, isotropic, multi-component Eddington inversion. It sees only the N-body
components themselves: **a disc is neither sampled nor included in its potential**, so the halo is in equilibrium only in
isolation and is perturbed when a disc is added (a true disc + halo equilibrium needs a made-to-measure code such as GalIC).
The gas disc is a perturbation on it, not a self-consistent partner.

**Gas halo** (`halo`). Komatsu & Seljak (2001) polytropic gas, hydrostatic in the dark halo, with Γ matched at r200. With the
default `consistent` mode the gas halo sits in the **same dark halo the disc and N-body parts use**; `--legacy` reproduces the
Fortran (which used a different scale length and mass normalisation). Two options change the halo so it does not overwhelm
the disc:
* `--halo-core F` adds a constant entropy floor K(r) = K_KS(r) + F·K_KS(r200) (K = P/ρ^{5/3}) and re-solves hydrostatic
  equilibrium keeping the gas mass inside r200 (F = 0 is exactly K&S). F = 0.1 raises the fraction of the default disc surviving
  pressure confinement from 19 % to 61 %; beyond ~3 R_d some truncation of a 10⁴ K disc is unavoidable.
* `--zero-edge` (after Donnert 2013) makes the gas hydrostatic only out to R_max = `rtrunc`·r200 with **P(R_max) = 0**, so the
  density and temperature fall smoothly to zero rather than ending in a pressure cliff. The mass inside r200 is matched by an
  overall entropy scale; the sampled halo then holds ≈ 1.01–1.03 × f_b M200 instead of up to 1.7×.

| halo (default disc) | disc mass kept | disc outer radius |
|---|---|---|
| K&S, free edge | 19 % | 2.6 kpc |
| K&S, zero edge | 17 % | 2.6 kpc |
| core 0.1, free edge | 61 % | 5.7 kpc |
| core 0.1, zero edge | 46 % | 4.5 kpc |

The Fortran also gives each halo particle a Gaussian velocity with σ = √(2u) per component (about the sound speed, ~250 km/s);
this is kept by default, `--cold` removes it. The spin parameter `lambda_b` is not applied by the Fortran and is ignored here.

**Compounding the disc and the halo.** Built independently, the two gas components overlap and the cold disc ignores the pressure of the
hot halo. `galaxy`/`system` therefore apply a pressure-balance rule: the disc exists where its pressure ρ_d c_s² is at least the
halo's P_h(r). The disc is cut at that height (rings whose midplane pressure is below the halo's vanish) and halo particles
inside the disc region are removed. `--naive-merge` just concatenates, as `merge_comps` did. This is approximate: the halo ignores
the disc's gravity and the disc ignores the halo's weight.

**AREPO background grid** (`--background-grid`; port of AREPO's `ADDBACKGROUNDGRID`, `src/add_backgroundgrid`). An oct-tree is
built on the gas; a massless cell is placed at each leaf-node centre; every particle's mass, thermal energy and momentum are spread
onto the cells within its smoothing length with the cubic-spline kernel (conserved exactly; the smoothing length is the distance to
the `--nngb`-th neighbour, since AREPO recomputes it). The box is `--boxsize-factor` (default 5, clamped to AREPO's 2–15) times the
gas extent, all species are shifted to the box centre, `BoxSize` is written to the header, and the file is double precision with a
`Masses` block (zero-mass cells are expected). Printed AREPO values: `BoxSize`, `ReferenceGasPartMass`; AREPO's `MeanVolume` comes
from its internal `GridSize`, which is not reproduced.

## Merging systems

A **system** is one galaxy (gas disc + gas halo + N-body parts, already compounded). `merge` reads a TOML file
([`examples/merger.toml`](examples/merger.toml)) listing any number of systems, each either **built** from a parameter file or
**loaded** from an existing GADGET HDF5 file (including Fortran output), and writes one stitched file.

```toml
output = "merger_ic.hdf5"
seed = 7
[background_grid]            # one grid on the combined gas
enabled = true
[gas_overlap]                # default mode = "pressure"
mode = "pressure"            # or "none"
k = 32
[orbit]                      # positions/velocities of systems 0 and 1 from a Kepler orbit (x-y plane, CM frame)
systems = [0, 1]
r_init = 700.0
r_peri = 30.0
ecc = 1.0
[[system]]
name = "primary"
params = "galaxy_params.txt"          # build ...
[system.options]                      # kwargs of build_system (halo_core, zero_edge, cold_halo, ...)
halo_core = 0.1
[system.overrides]                    # any Params field
ndisc = 40000
[[system]]
name = "secondary"
file = "other.hdf5"                   # ... or load
euler = [0, 60, 0]                    # z-x-z Euler angles [deg]
position = [300, 0, 0]                # (unless an [orbit] sets them)
velocity = [-50, 0, 0]
type_map = { 1 = 2 }                  # keep this system's halo as PartType2
taper_mass = 0.99                     # old merge_comps options: taper_mass,
reduce_stellar_mass = { type = 2, fraction = 0.5 }   #   reduce_stellar_mass,
excise_types = [3]                                   #   excise_types
```

**Gas interfaces between systems.** Where the gas of two systems overlaps, `merge` applies a pairwise pressure-balance rule in the
spirit of Donnert (2013), who resamples the overlapping gas rather than adding it: each gas particle's pressure
P = (γ−1)ρu is estimated from its `k` nearest neighbours, and a particle is removed if any other system's gas would have a higher
pressure at that point, so the higher-pressure gas wins. It needs only particles, so it works for systems loaded from files; the
price is a noisy estimate and a sharp cut. With `mode = "none"` the gas is superposed and a warning is printed if the gas extents
overlap. Choose `r_init` larger than the sum of the gas extents (`rtrunc·r200` each) for the cleanest start.

## Use the modules directly

The package is a set of independent modules, so any piece can be reused on its own. Self-contained, runnable scripts are in
[`examples/`](examples/) (they need the package installed, or `PYTHONPATH=.` from the repository root); they use
`galaxy_components.presets.example_params()` so no parameter file is needed.

| script | shows |
|---|---|
| [`01_gas_disc.py`](examples/01_gas_disc.py) | potential → vertical structure → rotation curve → particles |
| [`02_gas_halo.py`](examples/02_gas_halo.py) | K&S halo vs entropy core vs zero-pressure edge |
| [`03_nbody_halo.py`](examples/03_nbody_halo.py) | spherical N-body halo and a check against the model mass profile |
| [`04_background_grid.py`](examples/04_background_grid.py) | AREPO background grid for **any** gas particles (`--input` an HDF5 file) |
| [`05_stitch_systems.py`](examples/05_stitch_systems.py) | two galaxies, Kepler orbit, interface resolution, stitching |
| [`06_custom_potential.py`](examples/06_custom_potential.py) | compose potentials; reuse the hydrostatic solver for your own Φ(R,z) |

```python
import numpy as np
from galaxy_components.presets import example_params
from galaxy_components.system import build_system
from galaxy_components import merge

a = build_system(example_params(ndisc=20000), "a", halo_core=0.1)                  # one galaxy in memory
b = build_system(example_params(m200=30.0, ndisc=8000), "b", halo_core=0.1)
x1, v1, x2, v2 = merge.orbit_two_body(a.meta["m200"], b.meta["m200"], r_init=700, r_peri=30)
placed, _ = merge.resolve_gas_overlap([merge.place(a, x1, v1), merge.place(b, x2, v2, euler=(0, 60, 0))])
species, box = merge.finalise(merge.stitch(placed), background_grid=True)           # one background grid
merge.write_merged("merger.hdf5", species, box)
```

| module | contents |
|---|---|
| `potentials` | `Hernquist`, `NFW`, `LogarithmicHalo`, `ExponentialDisc` with `phi(R,z)`, `vc2(R)` |
| `model` | `build_model(params)` – derived radii/masses and the background components |
| `gas_disc` | `build_gas_disc`, `solve_rings` (any Φ), `sample_particles` |
| `gas_halo` | `build_gas_halo`, `add_entropy_core(..., zero_edge=)`, `sample_gas_halo` |
| `nbody` | `build_nbody`, `sample_nbody` |
| `compound` | `truncate_disc`, `excise_halo` (disc ↔ halo pressure balance) |
| `backgroundgrid` | `add_background_grid` (works on any particle arrays) |
| `system` | `build_system`, `read_system` |
| `merge` | `place`, `orbit_two_body`, `resolve_gas_overlap`, `stitch`, `finalise`, `merge_from_config` |
| `io`, `params`, `presets` | HDF5 writer, namelist reader, ready-made parameter sets |

## Tests

```bash
pytest -q              # everything, ~20 s
pytest -q -m smoke     # end-to-end runs of every command and example, ~13 s
```
The physics tests check against analytic results: the exponential-disc potential against the Freeman closed form, the vertical solver
against the isothermal sech² sheet, hydrostatic balance of every halo variant, conservation in the background grid, the Kepler orbit's
energy and angular momentum, and the N-body distribution function against a fresh run of the Fortran `mk_nbody_comps` (agreement < 1 %).
The smoke tests run each command and example on tiny inputs and check the files they write (consistent header, finite values,
unique IDs, particles inside the box). There is no Fortran oracle for the disc: the Fortran `mk_gas_disc` crashes with an
out-of-bounds read in `fintrp` (`legacy/mk_gas_disc.f90`, `xp(0)`).

## Differences from the Fortran

| area | Fortran | here |
|---|---|---|
| π | 3.14159 | exact |
| compile-time choices (`-DNFW_HALO`, …) | makefile flags | runtime (`&MODEL`) |
| gas halo potential | `a = r200/c200`, mass `M200` | same dark halo as the disc (`--legacy` restores) |
| N-body potential interpolation | extrapolates from the wrong table segment | interpolates correctly |
| N-body velocity scale | `√(G·mass_comp(1)/r_scale)` | `√(G·mtot/r_scale)` (identical for NFW/Burkert) |
| N-body `drptab` edge fill | mixes in `dptab` differences | same by default; `--fix-legacy` corrects |
| disc vertical positions | snapped to the 500-point z table | interpolated (`snap_z` restores) |
| disc rotation curve | last radius reads past the array (`i.lt.nr` bug) | fixed |
| bulge scale in rotation curve | `fbulge·rs`, inconsistent with its potential | same scale as the potential |
| halo recentring | shifts by the inner half's mean | not applied |
| output | Gadget binary or HDF5 | HDF5 only |

## Limitations

* No disc in the N-body potential (see above); systems are superposed, not re-equilibrated, apart from the gas pressure-balance rules.
* Glass-file sampling (`ispoisson=.false.`) and Gadget binary output are not ported; `form_comp='general'` is not supported.
* `read_tabulated` for a tabulated halo profile is a verbatim port that could not be tested (no data file).
* The AREPO output has not been run in AREPO.
* A 10⁴ K exponential disc inside a hot halo is always pressure-confined at a few scale lengths; see the table above.

## References

Komatsu & Seljak 2001, MNRAS 327, 1353 (gas halo) · Donnert 2013, arXiv:1311.7066 (finite-boundary gas, interface handling) ·
Hernquist 1990 · Springel, Di Matteo & Hernquist 2005 · Yurin & Springel 2014, MNRAS 444, 62 (GalIC) · Freeman 1970 (exponential
disc) · Weinberger, Springel & Pakmor 2020 (AREPO).
