"""Convert gas particles to an AREPO initial condition with a background grid.

Python port of AREPO's ``ADDBACKGROUNDGRID`` mode (src/add_backgroundgrid/):

  1. Build an oct-tree on the gas particles only.  A *node* exists wherever
     >= 2 particles share an octant of their parent (the root always exists);
     a *leaf* is a node none of whose children is a node.
  2. Put a new, massless gas cell at the centre of every leaf node (volume =
     node side^3).  Leaves are small where the gas is dense and large where
     it is sparse, so the cells tile the whole box adaptively.
  3. Spread each old particle's mass, thermal energy and momentum onto the
     new cells inside its smoothing length h with the cubic-spline kernel,
     weighted by cell volume and normalised so that each particle's
     contribution sums to exactly one (mass is conserved).
  4. Drop the old particles; u = E/m and v = P/m for cells with m > 0.

The box is set as in ``prepare_domain_backgroundgrid``: BoxSize is forced into
[2, 15] x the gas extent, and everything is shifted so the box is centred on
the gas.  AREPO does not write h to ICs; it recomputes it, so here h is the
distance to the ``nngb``-th nearest gas neighbour.
"""
import numpy as np
from scipy.spatial import cKDTree

FACTOR_MIN_BOX = 2.0
FACTOR_MAX_BOX = 15.0
MAXLEVEL = 21
# AREPO 3D cubic spline (support h)
C1, C2, C5 = 2.546479089470, 15.278874536822, 5.092958178941


def kernel(u, h):
    u = np.asarray(u)
    w = np.where(u < 0.5, C1 + C2 * (u - 1) * u * u, C5 * (1 - u) ** 3)
    return np.where(u < 1.0, w, 0.0) / h ** 3


def _spread(v):  # 21 bits -> every third bit
    v = v.astype(np.uint64) & np.uint64(0x1FFFFF)
    v = (v | (v << np.uint64(32))) & np.uint64(0x1F00000000FFFF)
    v = (v | (v << np.uint64(16))) & np.uint64(0x1F0000FF0000FF)
    v = (v | (v << np.uint64(8))) & np.uint64(0x100F00F00F00F00F)
    v = (v | (v << np.uint64(4))) & np.uint64(0x10C30C30C30C30C3)
    v = (v | (v << np.uint64(2))) & np.uint64(0x1249249249249249)
    return v


def morton_keys(pos, boxsize):
    n = 1 << MAXLEVEL
    ijk = np.clip((pos / boxsize * n).astype(np.int64), 0, n - 1)
    return _spread(ijk[:, 0]) << np.uint64(2) | _spread(ijk[:, 1]) << np.uint64(1) | _spread(ijk[:, 2])


def _unspread(v):
    v = v & np.uint64(0x1249249249249249)
    v = (v ^ (v >> np.uint64(2))) & np.uint64(0x10C30C30C30C30C3)
    v = (v ^ (v >> np.uint64(4))) & np.uint64(0x100F00F00F00F00F)
    v = (v ^ (v >> np.uint64(8))) & np.uint64(0x1F0000FF0000FF)
    v = (v ^ (v >> np.uint64(16))) & np.uint64(0x1F00000000FFFF)
    v = (v ^ (v >> np.uint64(32))) & np.uint64(0x1FFFFF)
    return v


def tree_leaves(pos, boxsize):
    """Centres and side lengths of the leaf nodes of the gas oct-tree."""
    keys = morton_keys(pos, boxsize)
    shift = lambda L: np.uint64(3 * (MAXLEVEL - L))
    nodes = {0: np.array([0], dtype=np.uint64)}
    for L in range(1, MAXLEVEL + 1):
        cells, counts = np.unique(keys >> shift(L), return_counts=True)
        nodes[L] = cells[counts >= 2]
    leaves = []
    for L in range(0, MAXLEVEL + 1):
        if L == MAXLEVEL:
            is_leaf = np.ones(nodes[L].size, dtype=bool)
        else:
            has_child = np.isin(nodes[L], np.unique(nodes[L + 1] >> np.uint64(3)))
            is_leaf = ~has_child
        ids = nodes[L][is_leaf]
        if ids.size:
            ix = _unspread(ids >> np.uint64(2)).astype(float)
            iy = _unspread(ids >> np.uint64(1)).astype(float)
            iz = _unspread(ids).astype(float)
            side = boxsize / (1 << L)
            leaves.append((np.column_stack((ix, iy, iz)) + 0.5) * side)
            leaves.append(np.full(ids.size, side))
    centres = np.concatenate(leaves[0::2]) if leaves else np.empty((0, 3))
    sides = np.concatenate(leaves[1::2]) if leaves else np.empty(0)
    return centres, sides


def smoothing_lengths(pos, nngb=64):
    k = min(nngb, len(pos))
    d, _ = cKDTree(pos).query(pos, k=k, workers=-1)
    return d[:, -1]


def _pairs(cell_tree, pos, h, chunk=20000):
    """Yield (i, j, r) for all cells j within h_i of particle i."""
    for s in range(0, len(pos), chunk):
        lists = cell_tree.query_ball_point(pos[s:s + chunk], r=h[s:s + chunk], workers=-1)
        counts = np.fromiter((len(l) for l in lists), dtype=np.int64, count=len(lists))
        i = np.repeat(np.arange(s, s + len(lists)), counts)
        j = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists]) if counts.sum() else np.empty(0, np.int64)
        yield i, j


def add_background_grid(pos, vel, mass, u, boxsize_factor=5.0, boxsize=None, nngb=64,
                        hsml=None, verbose=True):
    """Gas particles -> AREPO cells.

    pos, vel (N,3), mass (N,), u (N,) specific internal energy.  Returns a dict:
      pos, vel, mass, u, volume (new cells), offset (add to other species),
      boxsize, and AREPO parameter suggestions.
    """
    log = print if verbose else (lambda *a, **k: None)
    pos = np.asarray(pos, dtype=np.float64)
    vel = np.asarray(vel, dtype=np.float64)
    mass = np.broadcast_to(np.asarray(mass, dtype=np.float64), (len(pos),)).copy()
    u = np.broadcast_to(np.asarray(u, dtype=np.float64), (len(pos),)).copy()

    lo, hi = pos.min(0), pos.max(0)
    len_gas = float((hi - lo).max())
    bs = boxsize if boxsize else boxsize_factor * len_gas
    bs = float(np.clip(bs, FACTOR_MIN_BOX * len_gas, FACTOR_MAX_BOX * len_gas))
    offset = 0.5 * bs - 0.5 * (lo + hi)
    p = pos + offset
    log(f"gas extent {len_gas:.4g}, BoxSize {bs:.4g} (allowed {FACTOR_MIN_BOX}-{FACTOR_MAX_BOX} x extent)")

    centres, sides = tree_leaves(p, bs)
    vol = sides ** 3
    inside = np.all((centres > 0) & (centres < bs), axis=1)
    centres, sides, vol = centres[inside], sides[inside], vol[inside]
    log(f"{len(centres)} leaf cells for {len(p)} gas particles; side range "
        f"{sides.min():.3g} .. {sides.max():.3g}")

    h = smoothing_lengths(p, nngb) if hsml is None else np.asarray(hsml, dtype=float)
    ctree = cKDTree(centres)

    # pass 1: normalisation W_i = sum_j V_j w(r_ij/h_i)
    W = np.zeros(len(p))
    orphan_nn = None
    for i, j in _pairs(ctree, p, h):
        r = np.linalg.norm(p[i] - centres[j], axis=1)
        np.add.at(W, i, vol[j] * kernel(r / h[i], h[i]))
    orphan = W <= 0
    if orphan.any():
        log(f"warning: {orphan.sum()} particles have no cell within h; assigned to the nearest cell")
        orphan_nn = ctree.query(p[orphan])[1]

    # pass 2: distribute
    m_c = np.zeros(len(centres))
    e_c = np.zeros(len(centres))
    mom_c = np.zeros((len(centres), 3))
    for i, j in _pairs(ctree, p, h):
        keep = ~orphan[i]
        i, j = i[keep], j[keep]
        r = np.linalg.norm(p[i] - centres[j], axis=1)
        wgt = vol[j] * kernel(r / h[i], h[i]) / W[i]
        np.add.at(m_c, j, mass[i] * wgt)
        np.add.at(e_c, j, mass[i] * u[i] * wgt)
        for k in range(3):
            np.add.at(mom_c[:, k], j, mass[i] * vel[i, k] * wgt)
    if orphan.any():
        idx = np.where(orphan)[0]
        np.add.at(m_c, orphan_nn, mass[idx])
        np.add.at(e_c, orphan_nn, mass[idx] * u[idx])
        for k in range(3):
            np.add.at(mom_c[:, k], orphan_nn, mass[idx] * vel[idx, k])

    pos_ok = m_c > 0
    u_c = np.zeros(len(centres))
    v_c = np.zeros((len(centres), 3))
    u_c[pos_ok] = e_c[pos_ok] / m_c[pos_ok]
    v_c[pos_ok] = mom_c[pos_ok] / m_c[pos_ok][:, None]
    log(f"mass before {mass.sum():.8g}, after {m_c.sum():.8g}; {int((~pos_ok).sum())} empty (zero-mass) cells")

    return {
        "pos": centres, "vel": v_c, "mass": m_c, "u": u_c, "volume": vol,
        "offset": offset, "boxsize": bs, "hsml": h,
        "suggest": {"BoxSize": bs, "ReferenceGasPartMass": float(mass.sum() / len(mass)),
                    "median_leaf_volume": float(np.median(vol))},
    }
