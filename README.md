# build_galaxy_components

This repository contains a set of four codes that have been developed by [A/Prof Chris Power (ICRAR/UWA)](https://www.icrar.org/people/cpower/). The purpose of each is highlighted below:

* [`mk_gas_disc.exe`](mk_gas_disc.f90) - this code is used to compute the hydrostatic equilibrium gas distribution of an isothermal disc embedded within a gravitational potential, usually a composite of a stellar bulge, stellar disc, and dark matter halo. 

* [`mk_gas_halo.exe`](mk_gas_halo.f90) - this code is used to compute the hydrostatic equilibrium gas distribution in a gaseous halo, based on the formalism set out in Komatsu & Seljak 2001.

* [`mk_nbody_comps.exe`](mk_nbody_comps.f90) - this code is used to generate equilibrium n-body components of galaxies with arbitary numbers of components. (...)

* [`merge_comps.exe`](merge_comps.f90) - this code is used to merge gas components with N-body components.

> **Layout.** The original Fortran codes (sources and `makefile`) now live in [`legacy/`](legacy/); build and run
> them from there, e.g. `cd legacy && make`, then `./mk_gas_disc.exe ../parameters.txt`. The Python package
> [`galaxy_components/`](galaxy_components/) is the maintained implementation (see "Python port" below).
> `parameters.txt` stays at the top level because both read it.

## Usage

### Compilation

To use these codes, they first need to be compiled. Begin by modifying the makefile included in the repository. 
```
> vim makefile
```
There are a number of available options in the makefile that can be activated based on the user's requirements. The first of these catagorise the type of halo potential within which you would like to generate your galaxy componets:
```
#OPT += -DNFW_HALO                  % adding the potential of an analytic NFW halo
#OPT += -DHERNQUIST_HALO            % adding the potential of a Hernquist halo 
#OPT += -DLOGARITHMIC_HALO          % adding the potential of a logarithmic halo 
```
Only one of these halo options can be activiated at a time. 
```
#OPT += -DHERNQUIST_BULGE           % adding the potential of a Hernquist bulge
#OPT += -DSTELLAR_DISC_EXPONENTIAL  % adding the potential of an exponential stellar disk
```
These flags specify the potentials of other components in the N-body model that can optionally be included. Neither, either or both of these options can be activated at a time.
```
#OPT += -DEXPONENTIAL_DISC          % adding the potential of an exponential gas disk
```
This specifies the potential for the gas disc to be placed within the N-body components. This flag needs to be activated when running `mk_gas_disc.exe`.
```
#OPT += -DGalIC                     % turn on GalIC constants
```
This flag needs to be activated in order to use constants and definitions that are consistent with the GalIC N-body initial conditions generation code [(Yurin & Springel, 2014)](https://ui.adsabs.harvard.edu/abs/2014MNRAS.444...62Y/abstract). For example, in the case in which you would like to merge a GalIC N-body simulation and gas components generated with this code. 
```
#OPT += -DHDF5                      % use HDF5
```
This flag needs to be activated to read and generate HDF5 outputs. If turned on, it is necessary to specify the paths of the include and lib folders of your HDF5 installation within the makefile. 

Once you have modified the makefile for your system, the code can be compiled and the executable generated:
```
> make
```
This will compile all executable programs. If you would like to modify the executable for alternative combinations of components, run 
```
> make clean
```
and then recompile the executables with the modified makefile. 


### Inputs
```
> vim parameters.txt 
```
To specify the details of the potential within which the components are built, first modify the input parameters. The parameter.txt file contains a series of namelists; the format is as follows: 
```
&HALO
m200=1.e0               % Virial mass, M200, in units of 1e10 Msol
c200=10                 % NFW concentration
/

&DISC
mdisc=0.035             % Disc mass as a fraction of M200
fdisc=0.2               % Disc scale length as a fraction of halo scale radius
mstar=0.9               % The fraction of the disc mass that is contained in the stellar component
sigma0=0.0              % ?
fhole=0.0               % ?
fdref=-0.0              % ?
/

&BULGE
mbulge=0.0              % Bulge mass as a fraction of M200
fbulge=0.0              % Bulge scale radius as a fraction of the halo scale radius
/

&GAS_DISC
temp=1.d4               % Temperature of the gas disc in K
ndisc=200000            % Number of particles in the gas disc  
/

&GAS_HALO
fbaryon=0.16            % Baryon fraction, OmegaB/Omega0
nhalo=20000             % Number of particles in the gas halo
lambda_b=0.0            % ?
v_radial=0.0            % ?
/

&BLACK_HOLE
mbh=.0e-5               % Black hole mass as a fraction of M200
/

&NBODY      
ncomp=0                 % ?
/

&OUTPUT
disc_file='gas_disc.gdt' % Output gas disc Gadget file
halo_file='gas_halo.gdt' % Output gas halo Gadget file
snap_format=1            % Format of output file (1 - Gadget binary snap1, 2 - Gadget binary snap2, 3 - HDF5)
/

&MISC
glass_file='./glass_32.gdt' % Input glass file
ispoisson=.false            % Sample from a Poisson distribution?
/
```
For the `&MISC` options, only one of these options should be specified. To generate your equilibrium gas component, either:
* Sample your particles from a Poisson distribution, using the ran3() rountine from NR in F77.
or
* By replicating a cube of glass and deforming it to generate the required mass distribution. A glass file is provided in this repository for this use. 

### Running the codes

Each code for generating components accepts the parameters from the input file parameters.txt and can be operated in a similar manner:
```
> ./mk_gas_disc.exe ./parameters.txt    % to build a gas disc within the potential specified 

> ./mk_gas_halo.exe ./parameters.txt    % to build a gas halo within the potential specified

> ./mk_nbody_comps.exe ./parameters.txt % to build an N-body model with components specified 
```
To merge an N-body model with a gas component:
```
> ./merge_comps.exe -in ./gas_component -in ./Nbody_components -out ./merged_output
```
It  is necessary to list the gas component file as the first in the list. If the input/output files are HDF5 format, this needs to be specified by:
```
> ./merge_comps.exe -in_hdf5 ./gas_component.hdf5 -in_hdf5 ./Nbody_components.hdf5 -out_hdf5 ./merged_output.hdf5
```
These flags can be modified for Gadget inputs (-in_snap1 or -in_snap2) or Swift inputs (-in_swift). Similarly, outputs can be written in any of these formats (-out_snap1, -out_snap2 or -out_swift).

## Notes

The mean molecular weight is assumed to 4/(1+3*0.76)~1.21 appropriate for a gas disc; 
for the hot halo, it's assumed to be ionized.

At the moment, I assume the N-body component is set up with GalIC - so potentials and 
parameters are calculated to be consistent with this.


## Python port (`galaxy_components/`)

`mk_gas_disc` and `mk_nbody_comps` have a Python port that reads the same `parameters.txt`
(requires numpy, scipy, h5py, f90nml; tests with `pytest`):
```
> python -m galaxy_components disc   parameters.txt --outdir out
> python -m galaxy_components nbody  parameters.txt --outdir out
```
The makefile `-D` flags become an optional `&MODEL` namelist group (ignored by the Fortran
programs): `halo_type='hernquist'|'nfw'|'logarithmic'`, `galic=.true.`. Output is GADGET HDF5 only;
glass-file sampling (`ispoisson=.false.`) is not ported yet.

### What the N-body path does and does not model

The N-body sampler (`nbody`) is the default and only N-body path. It builds **spherical, isotropic
halo and bulge components** by Eddington inversion, each in the potential of the *N-body
components alone* (`&NBODY`: `ncomp`, `form_comp`, `fmass_comp`, ...). A bulge is simply another
component, e.g. `form_comp(2)='hernquist'`.

* **The stellar disc is not sampled.** `nbody` produces no disc particles, and the disc's gravity is
  **not** included in the halo/bulge potential, so these components are in equilibrium only in
  isolation. Adding a disc perturbs them; expect adiabatic response/relaxation when evolved.
  Use it where that is acceptable (e.g. a small disc-to-halo mass ratio), or provide the stellar
  disc and a consistent halo from a made-to-measure code such as GalIC.
* **The gas disc is a test-particle-like add-on.** `disc` solves hydrostatic equilibrium for the gas
  in the analytic halo + bulge + stellar-disc potential (the gas's own gravity enters only through
  the local vertical term, and its radial pull through the Freeman formula). It does not feed back
  on the halo, and the halo is not re-equilibrated to the disc.
* True disc + halo equilibrium (made-to-measure or Schwarzschild-type) is out of scope here.

### Gas halo, compounding with the disc, and the AREPO background grid

```
> python -m galaxy_components halo   parameters.txt --outdir out   # gas halo only (mk_gas_halo)
> python -m galaxy_components galaxy parameters.txt --outdir out --background-grid   # everything -> AREPO IC
```

**The gas halo** is the Komatsu & Seljak (2001) polytropic hydrostatic solution in a *dark-matter-only*
halo (`&GAS_HALO`: `fbaryon`, `nhalo`, `rtrunc`, `form`). Closed forms replace the Fortran quadratures;
it reproduces the Fortran gamma, T0 and radial profiles. Two things differ by default:
* `consistent` (default on, `--legacy` off): the Fortran builds the Hernquist gas halo with
  `a = r200/c200` and mass `M200`, whereas the disc and N-body codes use the GalIC scale
  `a = (r200/c200) sqrt(2 f(c))`. The gas halo therefore sat in a different potential from the disc.
  The default now uses the same dark halo everywhere (for c=10 this changes gamma 1.185 -> 1.143 and
  T0 3.8e6 -> 1.6e6 K).
* The Fortran gives each halo particle a Gaussian velocity of sigma = sqrt(2u) per component (~ the
  sound speed, ~250 km/s) on top of its thermal energy. This is kept by default; `--cold` removes it.
  `lambda_b` (spin) is not applied in the Fortran (code commented out) and is ignored here too.

**How the halo and disc compound.** The two gas components are computed independently, each assuming
the other is absent, and `merge_comps` simply concatenates them. Where they overlap this double counts
mass and ignores the pressure the hot halo exerts on the cold disc. `galaxy` instead (default) applies
a pressure-balance criterion: the disc exists where `P_disc = rho_d cs^2 >= P_halo(r)`. The disc is cut
at that height (rings whose midplane pressure is below the halo's vanish) and halo particles inside the
disc region are removed; `--naive-merge` gives the old behaviour. This is only approximately an
equilibrium: the halo ignores the disc's gravity and the disc ignores the halo's weight. With the
default parameters the effect is large (a 10^4 K disc is overpressured by a cuspy 1.5x10^6 K halo beyond
~1 R_d), which is a real property of those parameters rather than an artefact -- see the printed
`disc_mass_after_truncation`.

**Background grid** (`--background-grid`) is a port of AREPO's `ADDBACKGROUNDGRID`
(`src/add_backgroundgrid`): an oct-tree is built on the gas, a massless cell is placed at each leaf-node
centre, and every particle's mass, energy and momentum are spread onto the cells within its smoothing
length with the cubic-spline kernel (conserved exactly; the smoothing length is recomputed as the
distance to the `--nngb`-th neighbour, since AREPO does not read one from the IC). The box is set to
`--boxsize-factor` (default 5, clamped to AREPO's 2-15) times the gas extent, all species are shifted to
the box centre, and `BoxSize` is written to the header. Cells are written as double precision with a
`Masses` block. The AREPO values `BoxSize` and `ReferenceGasPartMass` are printed; AREPO's `MeanVolume`
comes from its internal `GridSize`, which is not reproduced here.

#### A halo that doesn't crush the disc: `--halo-core`

The K&S halo is cuspy, so its central pressure is high. `--halo-core f` (or `halo_core=f` in a
`&MODEL` group) keeps the same dark potential and the same gas mass inside r200 but adds a constant
entropy floor, K(r) = K_KS(r) + f K_KS(r200) with K = P/rho^(5/3), and re-solves hydrostatic
equilibrium (`gas_halo.add_entropy_core`; f = 0 is exactly K&S). Effect on the default parameters
(disc of 0.40 x 10^10 Msol, R_d = 2.1 kpc, T_disc = 10^4 K):

| f    | disc mass kept | disc outer radius | halo T(1 kpc) | halo n(1 kpc) |
|------|----------------|-------------------|---------------|---------------|
| 0    | 19 %           | 2.6 kpc           | 1.5e6 K       | 0.91 cm^-3    |
| 0.03 | 46 %           | 4.5 kpc           | 2.6e6 K       | 0.10 cm^-3    |
| 0.1  | 61 %           | 5.7 kpc           | 3.3e6 K       | 0.028 cm^-3   |
| 0.3  | 72 %           | 6.6 kpc           | 4.2e6 K       | 0.008 cm^-3   |
| 1    | 79 %           | 7.1 kpc           | 6.2e6 K       | 0.002 cm^-3   |

The core helps a lot but saturates: beyond ~3 R_d the *disc's own* midplane pressure (an exponential,
10^4 K gas) falls to ~10^3 K cm^-3, below that of any plausible hot halo at those radii, so some
truncation of the outer disc is unavoidable for this disc/halo pairing. Other ways to keep more disc: a
lower `fbaryon`, a hotter or more extended disc (`temp`, `fdisc`), or accepting an outer disc that is
pressure-confined. Note that a core pushes gas outwards: for fixed M_gas(<r200) there is more mass
between r200 and `rtrunc`*r200 (at f = 0.1 the sampled halo is ~1.6 x heavier), so reduce `rtrunc` if
that matters.

### Building systems and stitching them together (replaces `merge_comps`)

The old workflow built gas disc, gas halo and N-body parts as separate files and glued them with
`merge_comps.exe`. The Python workflow builds each **system** (one galaxy: gas disc + gas halo +
N-body components, already compounded by pressure balance) and then stitches any number of systems into one file.

```
> python -m galaxy_components system  galaxy.txt --outdir out --name gal1.hdf5   # a single system
> python -m galaxy_components merge   examples/merger.toml --outdir out          # stitch several
```
`merge` reads a TOML file (see `examples/merger.toml`). Each `[[system]]` is either built from a
parameter file (`params`, with optional `[system.overrides]` of any `&...` parameter and
`[system.options]` such as `halo_core`) or loaded from an existing GADGET HDF5 file (`file`, including
files from the Fortran codes). Per system you can set `position`, `velocity`, `euler = [phi, theta, psi]`
(degrees, z-x-z), `type_map = {1 = 2}` (keep a system's halo as a separate PartType) and the old
`merge_comps` tidying options: `taper_mass`, `reduce_stellar_mass = {type, fraction}`, `excise_types`.
An `[orbit]` table (`systems`, `r_init`, `r_peri`, `ecc`) sets the positions and velocities of two systems
on a Kepler orbit in the centre-of-mass frame. `[background_grid]` builds **one** AREPO background grid on
the combined gas (so all systems share one tree and one box), and everything is shifted into the box.

**Gas interfaces between systems.** Where the gas of two systems overlaps, `merge` applies a pairwise
pressure-balance rule (default; `[gas_overlap] mode = "pressure"`, `k = 32`) in the spirit of Donnert (2013),
who resamples the overlapping gas rather than adding it. Each gas particle's pressure
P = (gamma-1) rho u is estimated from its k nearest neighbours; a particle is removed if any other system's
gas would have a higher pressure at that point. So at an interface the higher-pressure gas wins and the
other is cut away, rather than the two being superposed (which double counts mass and leaves the boundary
out of pressure balance). It needs only particles, so it works for systems loaded from files; the cost is a
noisy estimate and a sharp cut. Set `mode = "none"` for plain superposition (a warning is printed if the gas
extents overlap). In the 400 kpc example only 0.2% of the secondary's gas is removed.

Remaining caveats: systems are not re-equilibrated (each N-body halo is sampled in isolation; the cut
interface is in pressure balance only to the accuracy of the k-neighbour estimate), and an interface cut
this way is not smoothed, so expect a short transient. Choose `r_init` larger than the sum of the gas
extents (`rtrunc * r200` each) for the cleanest start.

#### Hydrostatic to a finite edge with zero pressure: `--zero-edge`

Following Donnert (2013), `--zero-edge` (or `halo_zero_edge=.true.` in `&MODEL`) solves the halo gas
hydrostatically only out to R_max = `rtrunc` x r200 and imposes P(R_max) = 0 there, so the density and
temperature fall smoothly to zero at the edge instead of the gas ending in a pressure cliff. It combines
with `--halo-core`. With the edge pressure fixed, the gas mass inside r200 is matched by an overall entropy
scale instead; the sampled halo then holds ~1.01-1.03 x f_baryon M200 instead of 1.06-1.7 x (the free-edge
solution leaves extra mass between r200 and `rtrunc`). Pressure falls steeply near the edge
(P ~ (R_max - r)^2.5), so the effective gas edge is at ~0.9 R_max. Effect on the default disc:

| halo                | disc mass kept | disc outer radius | particles for nhalo = 20000 |
|---------------------|----------------|-------------------|-----------------------------|
| K&S, free edge      | 19 %           | 2.6 kpc           | 21241                       |
| K&S, zero edge      | 17 %           | 2.6 kpc           | 20206                       |
| core 0.1, free edge | 61 %           | 5.7 kpc           | 33775                       |
| core 0.1, zero edge | 46 %           | 4.5 kpc           | 20588                       |

The zero edge does not help the disc by itself (the edge pressure is irrelevant at a few kpc); it removes
the heavy low-density envelope, and with a core it moves mass back inside r200, which raises the central
pressure a little (hence 46 % rather than 61 %).
