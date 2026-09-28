# Data

All data live outside the repository under `LYALENSER_DATA` (default `/data/LyaLenser`, see `lyalenser/paths.py`).
Everything is public; `scripts/fetch_data.sh` (deltas, quasar catalogue, ACT, Planck) and
`scripts/fetch_lowz_catalogues.sh` (LSS catalogues, BOSS) download it with resumable `wget`
(the DESI server returns HTTP 503 under parallel downloads; the scripts retry).

## Raw inputs (`$LYALENSER_DATA/raw/`, about 110 GB with the ACT release)

| Path | Content | Used by |
|---|---|---|
| `desi/lya-deltas/delta-lya-0-0/Delta/delta-*.fits.gz` | DR1 Lya-forest deltas VAC v1.0, region A (1040-1205 A rest frame, 1028 HEALPix files, 5.8 GB). HDUs LAMBDA, METADATA (LOS_ID, RA, DEC in radians, Z), DELTA_BLIND (the blinding only affects picca's own distance conversion; the arrays are the measured fluctuations), WEIGHT, CONT. | `lyalenser.desi_io` |
| `desi/lya-deltas/delta-lyb-0-0/Delta/` | region B (920-1020 A, 300 files, 1.3 GB). Chunked independently of region A: match forests on LOS_ID, never on file number. | `lyalenser.desi_io` |
| `desi/qso_iron/QSO_cat_iron_cumulative_v0.fits` | DR1 quasar catalogue (every observation of every target; the reader keeps one row per TARGETID). | `lyalenser.qso_io` |
| `desi/lss_v1.5/` | DR1 LSS clustering catalogues v1.5: `LRG`, `ELG_LOPnotqso`, `QSO`, `BGS_BRIGHT-21.5`, per cap (NGC/SGC) with randoms 0-1 and `nz.txt`; weights WEIGHT (completeness, redshift failure, imaging systematics). | `scripts/lowz_catalogues.py` |
| `boss/` | BOSS DR12v5 CMASSLOWZTOT galaxies and random0 per cap (columns RA, DEC, Z, WEIGHT_SYSTOT/CP/NOZ; randoms carry Z). | `scripts/lowz_catalogues.py` |
| `act/baseline/` | ACT DR6 lensing baseline: kappa alm (lmax 4000), analysis mask (nside 4096), N_L, README. The mask is applied squared on the kappa side of cross-spectra, as the release prescribes. `act/release/` has all 13 variants (21 GB), `act/baseline/simulations/` 400 simulation alms (61 GB, unused now). | `lyalenser.cmb_maps` |
| `planck/PR4_variations/` | Planck PR4 lensing (Carron+2022, "2018-like" MV): `PR42018like_klm_dat_MV.fits`, mean field, `mask.fits.gz` (nside 2048), Galactic coordinates (rotated to equatorial on load). | `lyalenser.cmb_maps` |

## Products

| Path | Content | Made by |
|---|---|---|
| `stageb/basis_dr1_ab_z196.h5`, `stageb/basis_cross_dr1_z196.h5` | Kaiser basis correlations through the DESI windows and the continuum projection (forest-forest, quasar-forest). | stage 1 |
| `lowz_v5/` | Fiducial templates (iteration 15; the convergence maps are those of `lowz_v4` plus the derivative maps `dkappa_slice_*_alm.fits`, `dkappa_combined_alm.fits`; `lowz_v5_l1300/` the lmax = 1300 variant). | stage 2 |
| `lowz_v4/` | Iteration-14 templates: `kappa_slice_<z1>_<z2>_alm.fits`, `kappa_combined_alm.fits`, masks and unit-bias maps at nside 512, `summary.json` (tracer table, bias fits, spectra, Wiener weights per coverage class). `lowz_v4_l1300/` is the lmax = 1300 variant, `lowz_v3` the model-covariance Wiener version. | stage 2 |
| `stageb/dr1_lowz_v8/` | Fiducial forest auto-correlation fits (iteration 15): `fits.h5`, `dr1_lowz.json`/`.md`, `auto_subslabs.json`; `sightlines.h5`, `xi.h5` and `catalogue.h5` are hard links to `dr1_lowz_v7d`'s. | stages 5 and 6 |
| `stageb/dr1_qso_v2/` | Fiducial quasar-forest fits (iteration 15); `xi_qf.h5` and `catalogue.h5` hard-linked from `dr1_qso_v1d`. | stage 5 |
| `audit_v5/` | The 100-draw random-template nulls of iteration 15 (`null_{auto,cross}_0.json`). | stage 7 |
| `stageb/dr1_lowz_v7d/` | Iteration-14 forest auto-correlation run: `sightlines.h5` (3.9 GB, the `SightlineSet`), `xi.h5` (pair-count cells, fitted tables, cached redshift-resolved cells), `catalogue.h5` (pair accumulators, 0.9 GB), `fits.h5` (amplitude fits), `dr1_lowz.json`/`.md`, `auto_subslabs.json`. | stages 4 and 6 |
| `stageb/dr1_qso_v1d/` | Iteration-14 quasar-forest run: `xi_qf.h5`, `catalogue.h5` (1.2 GB), `fits.h5`, `dr1_qso.json`/`.md`. | stage 5 |
| `stageb/dr1_lowz_<tag>/`, `stageb/dr1_qso_<tag>/` | Robustness variants (`l500`, `l1300`, `rp20`, `rp40`, `kmed`, `kfine`; `_ws` suffix: the cross halves run on the workstation). | stage 8 |
| `mocks/` | Stage-A mock campaigns (32 GB, archived; `legacy/`). | not used |

The JSON summaries of the fiducial and robustness products are committed in `results/`
(`results/README.md`); the paper quotes only numbers that are in those files.

## Machines

- Workstation: `/data/LyaLenser` (raw data complete, all products; not synced).
- NERSC Perlmutter: `/global/cfs/cdirs/m4895/users/anze/LyaLenser_data` (raw complete, `stageb/` products of the
  NERSC runs; the quasar catalogue is a symlink to the public DR1 tree). See `docs/computing.md`.
