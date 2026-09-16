"""Why the DR1 injection-expectation slope moved from 1.014 to 1.041 with the corrected table.

The injection expectation (`inject.injection_test(..., expectation=True)`) replaces delta_p delta_q by the table's
xi at the TRUE (undisplaced) separation while the kernel and the mean field use the displaced geometry: with a
self-consistent table, xi' = d xi / d r_perp, it must return a slope of one, and it is the bookkeeping test of
signs, cuts and the band basis.

The iteration-8 table is deliberately NOT self-consistent in that sense: xi carries the same-wavelength term
N(r_perp) inside the first radial bin, because the mean field has to subtract it, while xi' omits it, because an
instrumental correlation between pixels at the same observed wavelength is not displaced by lensing. The
injection machinery, however, displaces the whole table, N included, so it injects a response the kernel does not
model and the recovered slope is biased high by exactly that term's share of the kernel.

This runs the same injection twice on the real DR1 sightlines: once with the table as the estimator uses it, and
once with N subtracted from xi (the lensable table). If the explanation is right the second returns to ~1.01.

Usage: python injection_same_wavelength.py --run $LYALENSER_DATA/stageb/dr1_lowz_v2 [--lowz $LYALENSER_DATA/lowz_split]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp

HERE = Path(__file__).resolve().parent
CODE = HERE.parent
for p in (CODE, CODE / 'pipeline', HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
from paths import DATA
from campaign4 import campaign_config, read_xi
from templates import sphere_band_templates
from inject import injection_test
from mock import load_sightlines
from xi_model import XiTable
from xi_spline import Envelope, XiCorrection
from xi_fit import BASIS


def lensable_table(table, basis_path, params):
    """Copy of ``table`` with the same-wavelength term removed from xi (xi_rp already excludes it)."""
    b = {k: read_xi(basis_path, f'projected/{k}') for k in BASIS}
    ref = XiTable(b[BASIS[0]].r_perp, b[BASIS[0]].r_par,
                  sum(x * b[k].xi for x, k in zip([1., 2 * 1.4, 1.4 ** 2], BASIS)),
                  sum(x * b[k].xi_rp for x, k in zip([1., 2 * 1.4, 1.4 ** 2], BASIS)), {})
    corr = XiCorrection(Envelope(ref))
    c = np.asarray(params['coefficients'], float)
    sw_only = np.concatenate([np.zeros(corr.n_s), c[corr.n_s:]])
    RP, RZ = np.meshgrid(table.r_perp, table.r_par, indexing='ij')
    return XiTable(table.r_perp, table.r_par, table.xi - corr.design(RP, RZ) @ sw_only, table.xi_rp,
                   dict(table.meta or {}, note='same-wavelength term removed from xi'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', type=Path, default=DATA / 'stageb/dr1_lowz_v2')
    ap.add_argument('--lowz', type=Path, default=DATA / 'lowz_split')
    ap.add_argument('--basis', type=Path, default=DATA / 'stageb/basis_dr1_v2.h5')
    ap.add_argument('--nside', type=int, default=512)
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    cfg = campaign_config(1.)
    t0 = time.perf_counter()
    sl = load_sightlines(a.run / 'sightlines.h5')
    table = read_xi(a.run / 'xi.h5', 'xi')
    params = json.loads(json.dumps(table.meta))['fit']
    alm = hp.read_alm(str(a.lowz / 'kappa_combined_alm.fits'))
    templates, _ = sphere_band_templates(alm, sl.ra, sl.dec, nside=a.nside, source='combined')
    alpha = sum(t.alpha for t in templates[:3])
    out = {}
    for label, tab in (('as used (xi carries N)', table),
                       ('lensable table (N removed from xi)', lensable_table(table, a.basis, params))):
        e = injection_test(sl, tab, alpha, [-.5, -.25, .25, .5], cfg, templates=templates, expectation=True)
        out[label] = {'paired_slopes_by_amplitude': e['paired_slopes_by_amplitude'], 'paired_slope': e['paired_slope']}
        print(f'[{time.perf_counter()-t0:6.0f} s] {label:36s} slopes {e["paired_slopes_by_amplitude"]}', flush=True)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(out, indent=1, default=float) + '\n')
        print('wrote', a.out)


if __name__ == '__main__':
    main()
