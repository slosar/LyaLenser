"""Size of the linear source-distance expansion of the lensing efficiency (paper, Section on the estimator).

The templates are built for a source at chi_ref (the weighted mean pixel distance of the forest sample); a pixel at
chi is displaced by alpha_T (1 + g1 (chi - chi_ref)) with one g1 for all templates (lyalenser.config.kernel_product_g1).
This script records the common g1 and effective lens distance, the per-slice values that a slice-by-slice treatment
would use (same kernel-product weight, restricted to the slice), the weighted distribution of pixel distances of the
production forest sample, and the resulting fractional change of the deflection across the sample.

python scripts/source_distance_expansion.py --run $LYALENSER_DATA/stageb/dr1_lowz_v7d --summary results/dr1_lowz_v7d.json
  -> results/source_distance_expansion.json (about 20 s; reads chi, w and slab of sightlines.h5 in chunks)
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py
from scipy.integrate import trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lyalenser.cosmo import chi as chi_of_z, z_of_chi, linear_pk_interp
from lyalenser.lensing import kernel, Z_CMB
from lyalenser.config import kernel_product_g1

SLICES = ((0.1, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.1), (1.1, 1.6))
ELLS = (40.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0)


def kernel_product_weight(cref):
    """The lens-distance weight of kernel_product_g1: kl x kCMB kernels times non-linear P, averaged over ELLS."""
    ccmb = float(chi_of_z(Z_CMB))
    chis = np.linspace(1.0, cref * (1.0 - 1e-5), 1200)
    zs = z_of_chi(chis)
    pk = linear_pk_interp(zmax=6.0, kmax=200.0, nonlinear=True)
    weight = np.zeros_like(chis)
    for ell in ELLS:
        weight += (2 * ell + 1) * kernel(chis, cref) * kernel(chis, ccmb) * pk.P(zs, (ell + .5) / chis, grid=False) / chis**2
    return chis, weight


def g1_of(mean_lens_chi, cref):
    return mean_lens_chi / (cref**2 * (1 - mean_lens_chi / cref))


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True); ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("results/source_distance_expansion.json"))
    a = ap.parse_args()
    summary = json.loads(a.summary.read_text()); cref = float(summary["chi_ref"])
    g1_common, lens_common = kernel_product_g1(cref)
    assert abs(g1_common - summary["config"]["g1"]) < 1e-9 * g1_common, "g1 differs from the production run"
    chis, weight = kernel_product_weight(cref)
    pix = pixel_distances(a.run)
    zmin, zmax = float(summary["zmin"]), float(summary["zmax"])
    c_lo, c_hi = float(chi_of_z(zmin)), float(chi_of_z(zmax))
    sig = pix["weighted_rms"]; dmean = pix["weighted_mean"] - cref

    def entry(g, lens):
        return dict(g1_per_mpch=float(g), effective_lens_chi=float(lens), effective_lens_z=float(z_of_chi(lens)),
                    rms_fraction=float(g * sig), factor_at_zmin=float(1 + g * (c_lo - cref)),
                    factor_at_zmax=float(1 + g * (c_hi - cref)),
                    mismatch_rms_fraction=float((g - g1_common) * sig),
                    # <A>/A_true = sum w G^2 (1+g1 dm)(1+g dm) / sum w G^2 (1+g1 dm)^2 to second order in g dm, with the
                    # pair weighting approximated by the pixel-weighted moments of dm = chi - chi_ref
                    amplitude_bias_estimate=float((g - g1_common) * (dmean + g1_common * (sig**2 + dmean**2))
                                                  / (1 + 2 * g1_common * dmean + g1_common**2 * (sig**2 + dmean**2))))

    out = dict(chi_ref=cref, z_ref=float(z_of_chi(cref)), zmin=zmin, zmax=zmax, chi_zmin=c_lo, chi_zmax=c_hi,
               ells=list(ELLS), pixel_distance=pix, common=entry(g1_common, lens_common), slices={})
    for z1, z2 in SLICES:
        m = (chis >= chi_of_z(z1)) & (chis < chi_of_z(z2))
        lens = float(trapezoid(weight[m] * chis[m], chis[m]) / trapezoid(weight[m], chis[m]))
        out["slices"][f"{z1:g}-{z2:g}"] = entry(g1_of(lens, cref), lens)
    a.out.write_text(json.dumps(out, indent=1))
    print(f"chi_ref={cref:.1f} (z={out['z_ref']:.3f}); pixel chi: mean {pix['weighted_mean']:.1f}, rms {sig:.1f}, "
          f"5-95% {pix['percentiles'][5]:.0f}-{pix['percentiles'][95]:.0f}")
    print(f"common: lens z={out['common']['effective_lens_z']:.2f} g1={g1_common:.3e} rms {out['common']['rms_fraction']:.3f} "
          f"ends {out['common']['factor_at_zmin']:.3f}/{out['common']['factor_at_zmax']:.3f}")
    for k, v in out["slices"].items():
        print(f"{k}: lens z={v['effective_lens_z']:.2f} g1={v['g1_per_mpch']:.3e} rms {v['rms_fraction']:.3f} "
              f"ends {v['factor_at_zmin']:.3f}/{v['factor_at_zmax']:.3f} mismatch rms {v['mismatch_rms_fraction']:+.3f} A bias {v['amplitude_bias_estimate']:+.4f}")


if __name__ == "__main__":
    main()
