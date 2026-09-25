"""Single place for machine-dependent paths. Override with the LYALENSER_DATA environment variable
(default /data/LyaLenser). Layout below it: raw/{desi,act,planck}, mocks/."""
import os
from pathlib import Path

DATA = Path(os.environ.get("LYALENSER_DATA", "/data/LyaLenser"))
RAW = DATA / "raw"
MOCKS = DATA / "mocks"
ACT_MASK = RAW / "act/baseline/mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits"
ACT_NL = RAW / "act/baseline/N_L_kk_act_dr6_lensing_v1_baseline.txt"
DESI_DELTAS = RAW / "desi/lya-deltas/delta-lya-0-0/Delta"
DELTA_ARRAYS = RAW / "desi/delta_forest_arrays.npz"
