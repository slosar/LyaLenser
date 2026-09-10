"""First look at the DESI DR1 Lya deltas: integrity, forest counts, redshift coverage, per-forest noise.
Writes report/data_delta_summary.json and report/figures/data_delta_summary.pdf. Standalone; reads
/data/LyaLenser/raw/desi/lya-deltas/delta-lya-0-0/Delta/*.fits.gz."""
import glob, json, sys, time
import numpy as np, healpy as hp
from astropy.io import fits
sys.path.insert(0, '/home/anze/Dropbox/work/LyaLenser/code')
from cosmo import chi as chi_of_z, hubble, C_KMS

LYA = 1215.67
files = sorted(glob.glob('/data/LyaLenser/raw/desi/lya-deltas/delta-lya-0-0/Delta/delta-*.fits.gz'))
ra, dec, zq, snr, npix, zmin_f, zmax_f = [], [], [], [], [], [], []
sig_delta_med, wmean = [], []
t0 = time.time(); bad = []
for i, f in enumerate(files):
    try:
        with fits.open(f, memmap=False) as h:
            lam = h['LAMBDA'].data; md = h['METADATA'].data
            d = h['DELTA_BLIND'].data; w = h['WEIGHT'].data
            ok = np.isfinite(d) & (w > 0)
            n = ok.sum(axis=1)
            # per-forest wavelength coverage
            lo = np.array([lam[r].min() if r.any() else np.nan for r in ok]); hi = np.array([lam[r].max() if r.any() else np.nan for r in ok])
            ra.append(md['RA']); dec.append(md['DEC']); zq.append(md['Z']); snr.append(md['MEANSNR']); npix.append(n)
            zmin_f.append(lo / LYA - 1); zmax_f.append(hi / LYA - 1)
            # rms of delta per forest (signal + noise) and mean weight
            dd = np.where(ok, d, 0.0)
            sig_delta_med.append(np.sqrt((dd ** 2).sum(axis=1) / np.maximum(n, 1)))
            wmean.append(np.where(n > 0, np.where(ok, w, 0).sum(axis=1) / np.maximum(n, 1), np.nan))
    except Exception as e:
        bad.append((f, str(e)))
    if i % 200 == 0: print(i, f'{time.time()-t0:.0f}s', flush=True)
ra, dec, zq, snr, npix = map(np.concatenate, (ra, dec, zq, snr, npix))
zmin_f, zmax_f, sig_delta_med, wmean = map(np.concatenate, (zmin_f, zmax_f, sig_delta_med, wmean))
print('bad files:', bad)
N = len(zq); print('forests', N)
# footprint area from nside-64 pixels containing >=1 forest (approximate)
pix = hp.ang2pix(64, np.degrees(ra), np.degrees(dec), lonlat=True)
area = len(np.unique(pix)) * hp.nside2pixarea(64, degrees=True)
# effective sightline density vs z: number of forests covering z per deg^2
zg = np.linspace(1.9, 3.6, 35); neff = [((zmin_f <= z) & (zmax_f >= z)).sum() / area for z in zg]
# noise: sigma_delta per pixel ~ 1/MEANSNR (flux S/N per pixel); P_N = sigma^2 dchi_pix with dchi_pix at the forest's mean z
zmid = 0.5 * (zmin_f + zmax_f); dchi = 0.8 / LYA * C_KMS / hubble(zmid) * 0.6766 * (1 + zmid) / (1 + zmid)  # Mpc/h per 0.8 A: c*dlambda/lambda /H * h ... (1+z) cancels
dchi = 0.8 * C_KMS / (LYA * (1 + zmid)) * (1 + zmid) / hubble(zmid) * 0.6766
sig = 1.0 / snr; PN = sig ** 2 * dchi
good = np.isfinite(PN) & (snr > 0)
out = dict(n_files=len(files), n_bad=len(bad), n_forests=int(N), area_deg2_nside64=float(area),
           z_quasar_pct=np.percentile(zq, [5, 25, 50, 75, 95]).tolist(),
           zmin_forest_pct=np.percentile(zmin_f[np.isfinite(zmin_f)], [5, 50, 95]).tolist(), zmax_forest_pct=np.percentile(zmax_f[np.isfinite(zmax_f)], [5, 50, 95]).tolist(),
           npix_pct=np.percentile(npix, [5, 50, 95]).tolist(),
           meansnr_pct=np.percentile(snr, [5, 25, 50, 75, 95]).tolist(),
           PN_pct=np.percentile(PN[good], [5, 25, 50, 75, 95]).tolist(), PN_median=float(np.median(PN[good])),
           lnPN_std=float(np.std(np.log(PN[good]))), harmonic_mean_ratio_k0p1={},
           rms_delta_pct=np.percentile(sig_delta_med[np.isfinite(sig_delta_med)], [5, 50, 95]).tolist(),
           neff_z=dict(z=zg.tolist(), n_per_deg2=[float(x) for x in neff]))
# harmonic-mean noise factor relative to equal weights at the median, for P1D = 0.19 and 0.10 Mpc/h
for p1d in (0.10, 0.19):
    eq = 1.0 / (np.median(PN[good]) + p1d); hm = np.mean(1.0 / (PN[good] + p1d))
    out['harmonic_mean_ratio_k0p1'][str(p1d)] = float(hm / eq)
json.dump(out, open('/home/anze/Dropbox/work/LyaLenser/report/data_delta_summary.json', 'w'), indent=1)
np.savez('/data/LyaLenser/raw/desi/delta_forest_arrays.npz', ra=ra, dec=dec, zq=zq, snr=snr, npix=npix, zmin_f=zmin_f, zmax_f=zmax_f, PN=PN, rms_delta=sig_delta_med, wmean=wmean, area_deg2=area)
print(json.dumps({k: v for k, v in out.items() if k != 'neff_z'}, indent=1))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
ax[0].hist(zq, bins=60); ax[0].set_xlabel('quasar z'); ax[0].set_title(f'{N} forests, ~{area:.0f} deg$^2$')
ax[1].plot(zg, neff); ax[1].set_xlabel('z'); ax[1].set_ylabel('forests covering z per deg$^2$'); ax[1].axhline(22, ls=':', c='k')
ax[2].hist(np.log10(PN[good]), bins=60); ax[2].set_xlabel(r'$\log_{10} P_N$ [Mpc/h] ($\sigma_\delta^2\Delta\chi$, $\sigma_\delta=1/$MEANSNR)'); ax[2].axvline(np.log10(0.33), c='r', ls='--', label='forecast 0.33'); ax[2].legend()
fig.tight_layout(); fig.savefig('/home/anze/Dropbox/work/LyaLenser/report/figures/data_delta_summary.pdf'); print('saved')
