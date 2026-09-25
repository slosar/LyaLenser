"""How much do the five tracer slices talk to each other in the pair estimator?

Each slice is currently fitted on its own: an N_t x N_t response matrix over that slice's science bands, curl
partners and junk, with the other slices absent from the model. Strictly the fit should be joint, one
(N_t x N_slice) matrix, because the slice templates are not orthogonal on the pair basis: they share a footprint,
the band filters leak across the cut sky, and neighbouring shells share large-scale modes.

This measures the leakage directly. For every ordered pair of slices it computes the full cross-response block
F[s,s'] = sum_pairs [ mm d_s d_s' + (mc/2)(d_s s_s' + d_s' s_s) + mcc s_s s_s' ]  (the same form `amplitude`
uses), then solves the slice-s fit against a unit science amplitude in slice s'. The result is the leakage matrix
L[s,s'] = d A_hat_s / d A_s' : the diagonal is one by construction, and row sums away from one are the bias the
separate fits carry when every slice really does have amplitude one.

Nothing is fitted to data here; only the response geometry is used, so the answer is noise-free.

Usage: python slice_crosstalk.py --run $LYALENSER_DATA/stageb/dr1_lowz_v4 [--lowz $LYALENSER_DATA/lowz_split]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp

ROOT = Path(__file__).resolve().parents[1]
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
from lyalenser.config import production_config
from lyalenser.templates import sphere_band_templates
from lyalenser.amplitude import pair_scalars, n_science
from lyalenser.desi_io import load_sightlines


def load_catalogue(path, group='all'):
    import h5py
    from lyalenser.pairs import PairCatalogue
    with h5py.File(path) as f:
        g = f[group]
        return PairCatalogue(g['a'][()], g['b'][()], g['thx'][()], g['thy'][()], g['theta'][()],
                             g['accum'][()].astype(np.float64), g['npair'][()], dict(g.attrs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', type=Path, default=DATA / 'stageb/dr1_lowz_v4')
    ap.add_argument('--lowz', type=Path, default=DATA / 'lowz_split')
    ap.add_argument('--nside-alpha', type=int, default=1024)
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    cfg = production_config()
    t0 = time.perf_counter()
    sl = load_sightlines(a.run / 'sightlines.h5')
    cat = load_catalogue(a.run / 'catalogue.h5')
    summary = json.loads((a.lowz / 'summary.json').read_text())
    names = [f"slice_{s['zmin']:g}_{s['zmax']:g}" for s in summary['slices']]
    print(f'{len(cat.a)} sightline pairs, {len(names)} slices', flush=True)

    bundles = []
    for name in names:
        alm = hp.read_alm(str(a.lowz / f'kappa_{name.replace("slice_", "slice_")}_alm.fits'))
        b, _ = sphere_band_templates(alm, sl.ra, sl.dec, nside=a.nside_alpha, source=name)
        bundles.append(b)
        print(f'[{time.perf_counter()-t0:6.0f} s] {name}: {len(b)} templates', flush=True)
    nt = len(bundles[0])
    nsci = sum(1 for t in bundles[0] if getattr(t, 'kind', '') == 'signal')

    x = cat.accum.sum(axis=2)
    g1 = cfg.g1
    mm = x[:, 3] + 2 * g1 * x[:, 4] + g1 * g1 * x[:, 5]
    mc = g1 * x[:, 6]
    mcc = g1 * g1 * x[:, 7]
    d_all = []; s_all = []
    for b in bundles:
        d, s, _ = pair_scalars(cat, b)
        d_all.append(d); s_all.append(s)
    d_all = np.concatenate(d_all, axis=0); s_all = np.concatenate(s_all, axis=0)
    print(f'[{time.perf_counter()-t0:6.0f} s] pair scalars {d_all.shape}', flush=True)

    F = (d_all * mm) @ d_all.T
    F = F + .5 * ((d_all * mc) @ s_all.T + (s_all * mc) @ d_all.T) + (s_all * mcc) @ s_all.T
    print(f'[{time.perf_counter()-t0:6.0f} s] full response {F.shape}', flush=True)

    ns = len(names)
    blocks = lambda i, j: F[i * nt:(i + 1) * nt, j * nt:(j + 1) * nt]
    # marginalised single-slice fit: one common science amplitude, every other component free
    M = np.zeros((nt, 1 + (nt - nsci)))
    M[:nsci, 0] = 1
    for j, i in enumerate(range(nsci, nt), 1):
        M[i, j] = 1
    e = np.zeros(nt); e[:nsci] = 1          # a unit common science amplitude in the source slice
    L = np.zeros((ns, ns)); rho = np.zeros((ns, ns))
    for i in range(ns):
        Fr = M.T @ blocks(i, i) @ M
        for j in range(ns):
            L[i, j] = np.linalg.solve(Fr, M.T @ (blocks(i, j) @ e))[0]
            rho[i, j] = (e @ blocks(i, j) @ e) / np.sqrt((e @ blocks(i, i) @ e) * (e @ blocks(j, j) @ e))
    np.set_printoptions(linewidth=200, precision=4, suppress=True)
    print('\ncommon-science response correlation between slices (rho):')
    print(rho)
    print('\nleakage matrix L[s,s\'] = dA_s/dA_s\' (rows: the fitted slice):')
    print(L)
    print('\nrow sums (what a separate per-slice fit returns when every slice truly has A = 1):')
    print(L.sum(axis=1))
    out = {'slices': names, 'rho': rho.tolist(), 'leakage': L.tolist(),
           'row_sums': L.sum(axis=1).tolist(), 'n_templates': int(nt), 'n_science': int(nsci),
           'note': 'L[s,s] = 1 by construction; row sums away from 1 are the bias of the separate per-slice fits'}
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(out, indent=1) + '\n')
        print('wrote', a.out)


if __name__ == '__main__':
    main()
