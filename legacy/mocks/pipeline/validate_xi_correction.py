"""Does the spline-corrected response table change the estimator's normalisation on mocks?

The iteration-8 table (`xi_spline`) fits ~60 extra parameters to the same pixel pairs that enter the estimator,
so it has two ways to go wrong that the two-parameter fit did not have: it can follow noise in the measured
correlation (a random error in the kernel, and a self-referential correlation between the kernel and the data it
multiplies), and it can distort the kernel where the data are weak. On the mocks the model shape is already
correct (the generator's own basis, the generator's own projection), so the correction has nothing real to fit:
any change in the paired response A(1) - A(0) is the cost of the extra freedom.

For each saved seed of an existing campaign this re-fits the response table both ways and refits the amplitudes
with the seed's own templates, and reports the paired response per method. It reads saved mocks only; it writes
nothing into the campaign products and performs no provenance check.

Usage:
  python validate_xi_correction.py --root $LYALENSER_DATA/mocks/iteration7/sparse --seeds 4008 4009 4010 \
      --basis $LYALENSER_DATA/mocks/iteration7/basis/basis.h5 --out ../../report/xi_correction_mocks.json
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
for p in (HERE, HERE.parent):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import run_mock_validation as v
from mock import sightlines_for_variant
from pairs import find_pairs
from xi_model import xi_from_data
from campaign4 import campaign_config, read_xi
from lowz import lowz_bundles


def load_basis(path):
    import h5py
    from xi_fit import BASIS
    out = {'projected': {k: read_xi(path, f'projected/{k}') for k in BASIS}}
    with h5py.File(path) as f:
        out['coarse'] = {k: f[f'coarse/{k}'][()] for k in BASIS}
    return out


VARIANTS = {'none':        dict(xi_correction='none'),
            'spline':      dict(xi_correction='spline'),
            'spline_nosw': dict(xi_correction='spline', xi_same_wavelength=False),
            'spline_medium':dict(xi_correction='spline', xi_knots='medium'),
            'spline_fine':dict(xi_correction='spline', xi_knots='fine'),
            'spline_small':dict(xi_correction='spline', xi_knots='small'),
            'spline_coarse':dict(xi_correction='spline', xi_knots='coarse')}


def one_seed(seed, root, basis, cfg, names=('combined', 'truth'), methods=('none', 'spline')):
    t0 = time.perf_counter()
    m = v.load_mock(Path(root) / f'seed{seed:d}.h5')
    b = lowz_bundles(m, cfg)
    pairs = find_pairs(m.sightlines, cfg.r_perp_max / m.sightlines.chi.min())
    out = {'seed': seed, 'methods': {}}
    counts = {}
    for A in (0, 1):
        sl = sightlines_for_variant(m, A, True)
        counts[A] = (sl, xi_from_data(sl, cfg).counts)
    for method in methods:
        c = cfg.copy(**VARIANTS[method])
        res = {}
        for A in (0, 1):
            sl, nd = counts[A]
            ft = v.table_for(sl, c, basis, counts=nd)
            cat = v.cat_for(sl, ft.table, c, pairs)
            from lowz import joint_amplitudes
            fits, jk, reg = joint_amplitudes(cat, b, c, sl, list(names))
            res[f'A{A}'] = {n: fits[n]['A'] for n in names}
            res.setdefault('xi_fit', {})[f'A{A}'] = {k: ft.params[k] for k in ('b_F2', 'beta_F', 'chi2') if k in ft.params}
        res['paired_response'] = {n: res['A1'][n] - res['A0'][n] for n in names}
        out['methods'][method] = res
        print(f'  seed {seed} {method:7s} paired response ' +
              ' '.join(f'{n}={res["paired_response"][n]:+.4f}' for n in names), flush=True)
    out['wall_s'] = time.perf_counter() - t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--basis', type=Path, required=True)
    ap.add_argument('--seeds', type=int, nargs='+', required=True)
    ap.add_argument('--out', type=Path, default=None)
    ap.add_argument('--scale', type=float, default=1.0)
    ap.add_argument('--methods', nargs='+', default=['none', 'spline'], choices=list(VARIANTS))
    a = ap.parse_args()
    cfg = campaign_config(a.scale)
    basis = load_basis(a.basis)
    results = []
    for s in a.seeds:
        d = Path(a.root) / str(s)
        try:
            results.append(one_seed(s, d if d.exists() else a.root, basis, cfg, methods=tuple(a.methods)))
        except Exception as exc:                                  # a missing seed must not kill the run
            print(f'  seed {s} FAILED: {exc}', flush=True)
    summary = {}
    for method in a.methods:
        for n in ('combined', 'truth'):
            vals = np.array([r['methods'][method]['paired_response'][n] for r in results if method in r['methods']])
            if len(vals):
                summary[f'{method}/{n}'] = {'n': len(vals), 'mean': float(vals.mean()),
                                            'sem': float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else None,
                                            'values': vals.tolist()}
    paired = {}
    for n in ('combined', 'truth'):
        for method in a.methods:
            if method == 'none':
                continue
            d = np.array([r['methods'][method]['paired_response'][n] - r['methods']['none']['paired_response'][n]
                          for r in results if method in r['methods'] and 'none' in r['methods']])
            if len(d):
                paired[f'{method}/{n}'] = {'mean_difference': float(d.mean()),
                                           'sem': float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else None}
    out = {'seeds': a.seeds, 'results': results, 'summary': summary, 'spline_minus_none': paired}
    print(json.dumps({'summary': summary, 'spline_minus_none': paired}, indent=1))
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(out, indent=1, default=float) + '\n')
        print('wrote', a.out)


if __name__ == '__main__':
    main()
