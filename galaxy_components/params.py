"""Read the Fortran-style namelist parameter file.

The existing ``parameters.txt`` is read unchanged.  The compile-time options of
the Fortran code (-DNFW_HALO, -DHERNQUIST_BULGE, -DGalIC, ...) become runtime
options in an optional extra ``&MODEL`` group, which the Fortran programs
ignore.  Defaults reproduce the makefile as shipped.
"""
from dataclasses import dataclass, field
from pathlib import Path

import f90nml

NCOMP_MAX = 5


@dataclass
class Params:
    # &HALO
    m200: float = 0.0
    v200: float = 0.0
    c200: float = 10.0
    sigma200: float = 0.0
    qflat: float = 1.0
    rcore: float = 0.0
    # &DISC
    mdisc: float = 0.0
    sigma0: float = 0.0
    fdisc: float = 1.0
    mstar: float = 0.0
    fhole: float = 0.0
    fdref: float = 0.0
    # &BULGE
    mbulge: float = 0.0
    fbulge: float = 0.0
    # &GAS_DISC
    temp: float = 1e4
    ndisc: int = 0
    # &BLACK_HOLE
    mbh: float = 0.0
    # &NBODY
    ncomp: int = 1
    npart_comp: list = field(default_factory=lambda: [0] * NCOMP_MAX)
    fmass_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)
    rscale_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)
    raniso_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)
    rtrunc_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)
    form_comp: list = field(default_factory=lambda: [""] * NCOMP_MAX)
    istaper: bool = False
    fdecay: float = 3000.0
    # &OUTPUT
    disc_file: str = "gas_disc"
    halo_file: str = "gas_halo"
    nbody_file: str = "nbody"
    snap_format: int = 3
    # &MISC
    glass_file: str = ""
    ispoisson: bool = True
    tabfile: str = ""
    istabulated: bool = False
    # &MODEL (Python only; replaces the Fortran -D flags)
    halo_type: str = "hernquist"  # 'hernquist' | 'nfw' | 'logarithmic'
    galic: bool = True  # GalIC conventions (total halo mass == mhalo)
    exponential_disc: bool = True
    # Gas-disc mesh (hardcoded in mk_gas_disc.f90)
    nrbin: int = 150
    nzbin: int = 500
    log_zmin: float = -15.0
    log_zmax: float = 1.3
    snap_z: bool = False  # snap vertical positions to the z table, as Fortran


_GROUPS = ("halo", "disc", "bulge", "gas_disc", "black_hole", "nbody", "output", "misc", "model")
_LIST_FIELDS = ("npart_comp", "fmass_comp", "rscale_comp", "raniso_comp", "rtrunc_comp", "form_comp")


def read_params(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"parameter file does not exist: {path}")
    nml = f90nml.read(str(path))
    p = Params()
    for group in _GROUPS:
        if group not in nml:
            continue
        for key, val in nml[group].items():
            key = key.lower()
            if not hasattr(p, key):
                continue
            if key in _LIST_FIELDS:
                cur = list(getattr(p, key))
                vals = val if isinstance(val, list) else [val]
                for i, v in enumerate(vals[:NCOMP_MAX]):
                    if v is not None:
                        cur[i] = v.strip().lower() if isinstance(v, str) else v
                setattr(p, key, cur)
            else:
                setattr(p, key, val.strip() if isinstance(val, str) else val)
    p.halo_type = p.halo_type.lower()
    return p


def output_name(base, snap_format):
    suffix = {1: ".snp1.gdt", 2: ".snp2.gdt", 3: ".gdt.hdf5", 4: ".swft.hdf5"}
    return base + suffix.get(snap_format if snap_format <= 4 else 1, "")
