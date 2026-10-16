"""Command-line entry points:  python -m galaxy_components {disc,nbody} parameters.txt"""
import argparse
from pathlib import Path

import numpy as np

from . import backgroundgrid, compound, gas_disc, gas_halo, nbody
from .constants import PI, GGDT
from .io import write_hdf5
from .model import build_model
from .params import output_name, read_params


def run_disc(args):
    p = read_params(args.params)
    model = build_model(p)
    disc = gas_disc.build_gas_disc(p, model)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    _write_disc_tables(out, p, model, disc)
    rng = np.random.default_rng(args.seed)
    pos, vel, ug, mpart = gas_disc.sample_particles(p, disc, p.ndisc, rng)
    species = {0: dict(pos=pos, vel=vel, ug=ug, mass=mpart)}
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)
    write_hdf5(out / output_name(p.disc_file, p.snap_format), species)


def _write_disc_tables(out, p, model, disc):
    r = disc.rtab
    hdr = (f"rdisc [kpc]: {model.rdisc:.8f}\nmdisc [1e10 Msol]: {model.mdisc:.8f}\n"
           f"mstar [1e10 Msol]: {model.mstar:.8f}\nsigma0 [1e10 Msol/kpc^2]: {disc.sigma0:.8f}\n"
           "Radius[kpc] Vc_bg[km/s] Vdisc[km/s] Vp[km/s] Vrot[km/s]")
    v = disc.vcomp
    np.savetxt(out / "vrot_components.txt",
               np.column_stack((r, np.sqrt(v["background"]), np.sqrt(v["gas"]),
                                np.sqrt(np.abs(v["pressure"])), disc.vrot)), header=hdr)
    dv = np.gradient(disc.vrot, r)
    kappa = np.sqrt(np.maximum(2 * (disc.vrot / r) * (dv + disc.vrot / r), 0.0))
    sig = disc.sigma0 * np.exp(-r / model.rdisc)
    q = disc.cs * kappa / PI / GGDT / sig
    np.savetxt(out / "disc_parameters.txt",
               np.column_stack((r, sig, disc.vrot, q, sig / 2 / disc.rho0)),
               header=hdr.split("\n")[0] + "\nRadius Sigma(R) Vrot Q Sigma/2rho0")


def run_nbody(args):
    p = read_params(args.params)
    nb = nbody.build_nbody(p, legacy=not args.fix_legacy)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rscale0 = nb.rscale_comp[0]
    np.savetxt(out / "distribution_function.txt",
               np.column_stack((nb.etab, rscale0 ** 3 * nb.ftab.T)),
               header="Binding Energy E [g*mtot/rscale]  F(E)")
    rng = np.random.default_rng(args.seed)
    comps = nbody.sample_nbody(p, nb, rng)
    s = nbody.summarise(nb, comps)
    print(f"Generated halo with {sum(len(c['pos']) for c in comps)} particles; {s['nvir']} inside r200")
    print("Centre of mass [kpc]:", s["com"], "\nCentre of mass velocity [km/s]:", s["vcom"])
    print("Velocity dispersion [km/s]:", s["sigma"], s["sigma_tot"])
    print("Radial velocity (dispersion) [km/s]:", s["vr"], s["sigma_r"])
    species = {i + 1: dict(pos=c["pos"], vel=c["vel"], mass=c["mpart"]) for i, c in enumerate(comps)}
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)
    write_hdf5(out / output_name(p.nbody_file, p.snap_format), species)


def run_halo(args):
    p = read_params(args.params)
    halo = gas_halo.build_gas_halo(p, consistent=not args.legacy)
    rng = np.random.default_rng(args.seed)
    pos, vel, ug, mpart = gas_halo.sample_gas_halo(p, halo, rng, velocity_noise=not args.cold)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    species = {0: dict(pos=pos, vel=vel, ug=ug, mass=mpart)}
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)
    write_hdf5(out / output_name(p.halo_file, p.snap_format), species)


def run_galaxy(args):
    """Disc + halo gas (+ N-body) in one file, optionally with an AREPO background grid."""
    p = read_params(args.params)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    model = build_model(p)
    disc = gas_disc.build_gas_disc(p, model, verbose=False)
    halo = gas_halo.build_gas_halo(p, consistent=not args.legacy_halo)
    hpos, hvel, hug, hmp = gas_halo.sample_gas_halo(p, halo, rng, velocity_noise=not args.cold_halo)

    if args.naive_merge:
        sdisc = disc
    else:
        sdisc, mfrac, _ = compound.truncate_disc(disc, halo)
        print(compound.mass_summary(disc, sdisc, mfrac))
        keep = compound.excise_halo(hpos, disc, halo)
        print(f"excised {int((~keep).sum())} of {len(keep)} halo particles inside the disc")
        hpos, hvel, hug = hpos[keep], hvel[keep], hug[keep]
    dpos, dvel, dug, dmp = gas_disc.sample_particles(p, sdisc, p.ndisc, rng)

    gas = dict(pos=np.concatenate((dpos, hpos)), vel=np.concatenate((dvel, hvel)),
               ug=np.concatenate((dug, hug)),
               mass=np.concatenate((np.full(len(dpos), dmp), np.full(len(hpos), hmp))))
    species = {0: gas}

    if not args.no_nbody and p.ncomp > 0 and p.npart_comp[0] > 0:
        nb = nbody.build_nbody(p, verbose=False)
        for i, c in enumerate(nbody.sample_nbody(p, nb, rng, verbose=False)):
            species[i + 1] = dict(pos=c["pos"], vel=c["vel"], mass=c["mpart"])
    if p.mbh > 0:
        species[5] = dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)), mass=p.mbh)

    boxsize = 0.0
    if args.background_grid:
        bgo = backgroundgrid.add_background_grid(gas["pos"], gas["vel"], gas["mass"], gas["ug"],
                                                 boxsize_factor=args.boxsize_factor, nngb=args.nngb)
        species[0] = dict(pos=bgo["pos"], vel=bgo["vel"], mass=bgo["mass"], ug=bgo["u"])
        for t, s in species.items():
            if t != 0:
                s["pos"] = np.asarray(s["pos"], dtype=float) + bgo["offset"]
        boxsize = bgo["boxsize"]
        print("AREPO suggestions:", bgo["suggest"])
    fname = out / (args.name or "galaxy_ic.hdf5")
    write_hdf5(fname, species, boxsize=boxsize, double=args.background_grid)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="galaxy_components")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("disc", run_disc), ("nbody", run_nbody), ("halo", run_halo), ("galaxy", run_galaxy)):
        sp = sub.add_parser(name)
        sp.add_argument("params")
        sp.add_argument("--outdir", default=".")
        sp.add_argument("--seed", type=int, default=19721)
        if name == "nbody":
            sp.add_argument("--fix-legacy", action="store_true",
                            help="correct Fortran edge-extrapolation quirk instead of reproducing it")
        if name == "halo":
            sp.add_argument("--legacy", action="store_true",
                            help="Fortran halo potential (a=r200/c200, M200) instead of the disc's dark halo")
            sp.add_argument("--cold", action="store_true",
                            help="no Gaussian velocity noise (Fortran adds sigma=sqrt(2u))")
        if name == "galaxy":
            sp.add_argument("--legacy-halo", action="store_true")
            sp.add_argument("--cold-halo", action="store_true")
            sp.add_argument("--naive-merge", action="store_true",
                            help="plain concatenation (merge_comps behaviour): no pressure-balance")
            sp.add_argument("--no-nbody", action="store_true")
            sp.add_argument("--background-grid", action="store_true", help="AREPO background grid")
            sp.add_argument("--boxsize-factor", type=float, default=5.0)
            sp.add_argument("--nngb", type=int, default=64)
            sp.add_argument("--name")
        sp.set_defaults(fn=fn)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
