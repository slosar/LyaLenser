"""Build the DR1 response basis: Hankel-transformed basis spectra pushed through the DESI pixel windows and the
per-forest continuum projection averaged over REAL forest pairs.

The three basis spectra P_lin(k) F(k) mu^{2i} (i = 0, 1, 2) are transformed with the line-of-sight pixel and
resolution windows (`xi_fit.basis_tables`, provider 'hankel'), then projected. Iteration 5 projected through a
single representative forest with uniform weights; for DR1 that forest was the whole slab (1295 pixels of
712 Mpc/h) while the median picca forest is 689 pixels (397 Mpc/h), and the resulting model over-predicts xi by
5-20 per cent over the kernel range. Here the projection is the pair-weight average over sampled real forest
pairs (`xi_fit.project_fine_sample`): the projector uses every pixel picca fitted the continuum to, the pair sums
use only the pixels the measurement keeps.

Usage:  python build_basis_dr1.py --out $LYALENSER_DATA/stageb/basis_dr1_v2.h5 [--files 40] [--pairs 400]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import h5py
from astropy.io import fits

HERE = Path(__file__).resolve().parent
CODE = HERE.parent
for p in (CODE, CODE / 'pipeline', HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
from paths import DATA
from cosmo import chi as chi_of_z
from campaign4 import campaign_config
from xi_fit import basis_tables, project_fine_sample, refine_rp_grid, coarse_bin, BASIS
from desi_io import delta_files

LYA = 1215.67


def forest_geometries(zmin, zmax, n_files=40, min_pixels=50, seed=7, verbose=True):
    """(chi, weight, in_range) for every forest of a random sample of delta files.

    ``chi`` and ``weight`` cover the whole picca forest (every pixel with WEIGHT > 0), because the continuum was
    fitted to all of them; ``in_range`` marks the pixels the measurement keeps.
    """
    rng = np.random.default_rng(seed)
    files = delta_files()
    sub = [files[i] for i in rng.choice(len(files), min(n_files, len(files)), replace=False)]
    out = []
    for i, f in enumerate(sub):
        with fits.open(f, memmap=False) as h:
            lam = h['LAMBDA'].data
            z = lam / LYA - 1
            chi_l = np.asarray(chi_of_z(z), np.float64)
            inz = (z >= zmin) & (z <= zmax)
            d = h['DELTA_BLIND'].data
            w = h['WEIGHT'].data
            ok = np.isfinite(d) & (w > 0)
            for j in range(ok.shape[0]):
                m = ok[j]
                if (m & inz).sum() < min_pixels:
                    continue
                out.append((chi_l[m], w[j][m].astype(np.float64), inz[m]))
        if verbose and i % 10 == 0:
            print(f'  {i}/{len(sub)} files, {len(out)} forests', flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=DATA / 'stageb/basis_dr1_v2.h5')
    ap.add_argument('--files', type=int, default=40, help='delta files sampled for the forest geometries')
    ap.add_argument('--pairs', type=int, default=400, help='forest pairs averaged in the projection')
    ap.add_argument('--rp-step', type=float, default=0.5, help='r_perp step of the projection, refined to xi_step')
    ap.add_argument('--los-pixel', type=float, default=0.55, help='line-of-sight pixel size, Mpc/h')
    ap.add_argument('--los-resolution', type=float, default=0.4, help='line-of-sight resolution sigma, Mpc/h')
    ap.add_argument('--zmin', type=float, default=2.1)
    ap.add_argument('--zmax', type=float, default=3.0)
    ap.add_argument('--seed', type=int, default=11, help='forest-pair sampling seed')
    ap.add_argument('--geometry-seed', type=int, default=7, help='delta-file sampling seed')
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    cfg = campaign_config(1.)
    t0 = time.perf_counter()

    forests = forest_geometries(a.zmin, a.zmax, a.files, seed=a.geometry_seed)
    npix = np.array([len(f[0]) for f in forests])
    nin = np.array([int(f[2].sum()) for f in forests])
    print(f'{len(forests)} forests; picca length median {np.median(npix):.0f} pixels, '
          f'kept {np.median(nin):.0f} ({time.perf_counter()-t0:.0f} s)', flush=True)

    raw = basis_tables(cfg, 'hankel', nk=cfg.analytic_nk, model=cfg.forest_model,
                       los_pixel=a.los_pixel, los_resolution=a.los_resolution)
    print(f'basis spectra transformed ({time.perf_counter()-t0:.0f} s)', flush=True)
    rp_grid = np.arange(0, cfg.xi_max + .5 * a.rp_step, a.rp_step)
    proj = project_fine_sample(raw, forests, cfg, n_pairs=a.pairs, seed=a.seed, rp_grid=rp_grid, verbose=True)
    proj = refine_rp_grid(proj, cfg)
    coarse = {k: coarse_bin(proj[k], cfg) for k in BASIS}
    print(f'projection done ({time.perf_counter()-t0:.0f} s)', flush=True)

    if a.out.exists():
        a.out.unlink()
    for k in BASIS:
        raw[k].save(a.out, f'raw/{k}')
        proj[k].save(a.out, f'projected/{k}')
    with h5py.File(a.out, 'a') as f:
        if 'coarse' in f:
            del f['coarse']
        for k in BASIS:
            f[f'coarse/{k}'] = coarse[k]
        f.attrs['los_pixel'] = a.los_pixel
        f.attrs['los_resolution'] = a.los_resolution
        f.attrs['provider'] = 'hankel'
        f.attrs['projection'] = 'pair-weighted over real DR1 forest pairs'
        f.attrs['n_forest_pairs'] = a.pairs
        f.attrs['n_files'] = a.files
        f.attrs['n_forests'] = len(forests)
        f.attrs['median_picca_pixels'] = float(np.median(npix))
        f.attrs['median_kept_pixels'] = float(np.median(nin))
        f.attrs['forest_model'] = cfg.forest_model
        f.attrs['seed'] = a.seed
    print(f'wrote {a.out} ({time.perf_counter()-t0:.0f} s)')


if __name__ == '__main__':
    main()
