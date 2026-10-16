"""GADGET-style HDF5 output (mode 3 of io.f90::write_hdf5)."""
import numpy as np
import h5py

NTYPE = 6


def write_hdf5(path, species, boxsize=0.0, verbose=True):
    """Write particles.

    species: dict ptype(0..5) -> dict(pos, vel, mass(float, per-particle mass),
             ug (optional internal energy)).  IDs run consecutively over types.
    Masses are stored in the header MassTable when every present type has a
    positive (uniform) mass, as in the Fortran writer.
    """
    npart = np.zeros(NTYPE, dtype=np.int32)
    masstab = np.zeros(NTYPE)
    for t, s in species.items():
        npart[t] = len(s["pos"])
        masstab[t] = s["mass"]
    with h5py.File(path, "w") as f:
        h = f.create_group("Header")
        h.attrs["NumPart_ThisFile"] = npart
        h.attrs["NumPart_Total"] = npart
        h.attrs["NumPart_Total_HighWord"] = np.zeros(NTYPE, dtype=np.int32)
        h.attrs["MassTable"] = masstab
        h.attrs["Time"] = 0.0
        h.attrs["BoxSize"] = float(boxsize)
        h.attrs["Redshift"] = 1.0  # sic: matches Fortran gadget_header default
        h.attrs["Omega0"] = 1.0
        h.attrs["OmegaLambda"] = 0.0
        h.attrs["HubbleParam"] = 1.0
        h.attrs["NumFilesPerSnapshot"] = np.int32(1)
        for name in ("Flag_Cooling", "Flag_DoublePrecision", "Flag_Entropy_ICs", "Flag_Feedback",
                     "Flag_Metals", "Flag_Sfr", "Flag_StellarAge"):
            h.attrs[name] = np.int32(0)
        first = 1
        for t in sorted(species):
            s = species[t]
            n = len(s["pos"])
            g = f.create_group(f"PartType{t}")
            g.create_dataset("Coordinates", data=np.asarray(s["pos"], dtype=np.float32))
            g.create_dataset("Velocities", data=np.asarray(s["vel"], dtype=np.float32))
            g.create_dataset("ParticleIDs", data=np.arange(first, first + n, dtype=np.int32))
            if "ug" in s:
                g.create_dataset("InternalEnergy", data=np.asarray(s["ug"], dtype=np.float32))
            first += n
    if verbose:
        print(f"Wrote {path}: N = {npart.tolist()}")
