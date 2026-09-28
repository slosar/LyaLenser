"""Size of the source-distance dependence of the templates (paper, Section on the templates).

The templates are built for a source at chi_ref (the weighted mean pixel distance of the forest sample); a pixel at
chi is displaced by alpha + (chi - chi_ref) dalpha, with dalpha the derivative map of the template
(lowz_catalogues.py, dW/dchi_s in place of W). This script records, per slice and for the combined template, the
effective scalar ratio r = sum_l (2l+1) C_l^{kappa dkappa} / sum_l (2l+1) C_l^{kappa kappa} over the science range
(the g1 a scalar treatment would use), the weighted distribution of pixel distances of the production forest sample,
and the resulting fractional change of the deflection across the sample (r sigma_chi rms; the factors at the ends of
the forest range). The common value used in iterations 11-14 (one CMB-kernel-weighted coefficient for every
template) is kept for reference.

python scripts/source_distance_expansion.py --run $LYALENSER_DATA/stageb/dr1_lowz_v8 --lowz $LYALENSER_DATA/lowz_v5 \
    --summary results/dr1_lowz_v8.json --previous results/dr1_lowz_v7d.json -> results/source_distance_expansion.json
(about 20 s; reads chi, w and slab of sightlines.h5 in chunks)
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py, healpy as hp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lyalenser.cosmo import chi as chi_of_z, z_of_chi

BAND_EDGES = (40, 200, 400, 600, 800, 1000)


def pixel_distances(run, chunk=50_000_000):
    """Weighted moments of the pixel distance over the production sightline set (all pixels with slab >= 0)."""
    sw = swc = swc2 = 0.0
    hist = None
    edges = np.linspace(3000., 5000., 2001)
    with h5py.File(Path(run) / "sightlines.h5", "r") as f:
        g = f["sightlines"]; n = g["chi"].shape[0]
        for i in range(0, n, chunk):
            c = g["chi"][i:i + chunk].astype(np.float64); w = g["w"][i:i + chunk].astype(np.float64)
            s = g["slab"][i:i + chunk] >= 0
            c, w = c[s], w[s]
            sw += w.sum(); swc += (w * c).sum(); swc2 += (w * c * c).sum()
            h, _ = np.histogram(c, bins=edges, weights=w); hist = h if hist is None else hist + h
    mean = swc / sw; rms = np.sqrt(swc2 / sw - mean**2)
    cdf = np.cumsum(hist) / hist.sum(); mid = .5 * (edges[1:] + edges[:-1])
    pct = {p: float(np.interp(p / 100., cdf, mid)) for p in (5, 16, 50, 84, 95)}
    return dict(weighted_mean=float(mean), weighted_rms=float(rms), n_pixels=int(n), percentiles=pct)


def ratios(lowz, name, edges=BAND_EDGES):
    """Effective derivative ratio of a template over the science range and per band, from the alms."""
    kk = hp.read_alm(str(Path(lowz) / f"kappa_{name}_alm.fits")); dk = hp.read_alm(str(Path(lowz) / f"dkappa_{name}_alm.fits"))
    ckk = hp.alm2cl(kk); ckd = hp.alm2cl(kk, dk); ell = np.arange(len(ckk)); w = 2 * ell + 1
    def r(lo, hi):
        m = (ell >= lo) & (ell < hi); return float(np.sum(w[m] * ckd[m]) / np.sum(w[m] * ckk[m]))
    return r(edges[0], edges[-1]), {f"{lo}-{hi}": r(lo, hi) for lo, hi in zip(edges[:-1], edges[1:])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True); ap.add_argument("--lowz", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True, help="dr1_lowz JSON of the run (chi_ref, derivative ratios at the sightlines)")
    ap.add_argument("--previous", type=Path, default=None, help="dr1_lowz JSON of an iteration 11-14 run (its common scalar g1, for reference)")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "results/source_distance_expansion.json")
    a = ap.parse_args()
    summary = json.loads(a.summary.read_text()); cref = float(summary["chi_ref"])
    pix = pixel_distances(a.run); sig = pix["weighted_rms"]
    zmin, zmax = float(summary["zmin"]), float(summary["zmax"]); c_lo, c_hi = float(chi_of_z(zmin)), float(chi_of_z(zmax))
    lowz = json.loads((a.lowz / "summary.json").read_text())
    names = ["combined"] + [f"slice_{s['zmin']:g}_{s['zmax']:g}" for s in lowz["slices"]]

    def entry(r):
        return dict(ratio_per_mpch=float(r), rms_fraction=float(r * sig), factor_at_zmin=float(1 + r * (c_lo - cref)), factor_at_zmax=float(1 + r * (c_hi - cref)))

    out = dict(chi_ref=cref, z_ref=float(z_of_chi(cref)), zmin=zmin, zmax=zmax, chi_zmin=c_lo, chi_zmax=c_hi, pixel_distance=pix,
               templates=a.lowz.name, science_range=[BAND_EDGES[0], BAND_EDGES[-1]], per_template={})
    for n in names:
        r, per_band = ratios(a.lowz, n)
        out["per_template"][n] = dict(entry(r), per_band={k: float(v) for k, v in per_band.items()},
                                      ratio_at_sightlines=float(summary.get("derivative_ratio", {}).get(n, float("nan"))))
    if a.previous is not None:
        prev = json.loads(a.previous.read_text()); g = float(prev["config"]["g1"])
        out["previous_common_scalar"] = dict(entry(g), source=str(a.previous), note="one CMB-kernel-weighted coefficient for all templates (iterations 11-14)")
    a.out.write_text(json.dumps(out, indent=1))
    print(f"chi_ref={cref:.1f} (z={out['z_ref']:.3f}); pixel chi: mean {pix['weighted_mean']:.1f}, rms {sig:.1f}, 5-95% {pix['percentiles'][5]:.0f}-{pix['percentiles'][95]:.0f}")
    for n, v in out["per_template"].items():
        print(f"{n:16s} ratio {v['ratio_per_mpch']:.3e} (at sightlines {v['ratio_at_sightlines']:.3e}) rms {v['rms_fraction']:.3f} ends {v['factor_at_zmin']:.3f}/{v['factor_at_zmax']:.3f}")
    if "previous_common_scalar" in out:
        v = out["previous_common_scalar"]; print(f"previous common  g1 {v['ratio_per_mpch']:.3e} rms {v['rms_fraction']:.3f} ends {v['factor_at_zmin']:.3f}/{v['factor_at_zmax']:.3f}")


if __name__ == "__main__":
    main()
