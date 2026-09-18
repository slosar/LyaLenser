"""CMB lensing convergence maps for the Stage B cross-checks: ACT DR6 baseline and Planck PR4 (Carron+2022,
'2018-like' MV), as HEALPix maps at the working nside with the mask each release prescribes for cross-spectra
(ACT: the analysis mask squared on the kappa side, per the release README; Planck: the release mask)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import healpy as hp
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline'):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import RAW

ACT=RAW/'act/baseline'
PLANCK=RAW/'planck/PR4_variations'
SURVEYS=('ACT','Planck')


MASKED_ON_INPUT={'ACT':True,'Planck':False}   # ACT README: treat the map as already carrying M^2; Planck: multiply by the release mask


def load_kappa(survey,nside,lmax):
    """(kappa_map, mask) of ``survey`` in ('ACT', 'Planck'); the alm are cut at ``lmax`` before the synthesis.
    Build the NaMaster field with ``masked_on_input=MASKED_ON_INPUT[survey]``."""
    if survey=='ACT':
        alm=np.nan_to_num(hp.read_alm(str(ACT/'kappa_alm_data_act_dr6_lensing_v1_baseline.fits')))   # NaN in unused modes
        mask=hp.ud_grade(hp.read_map(str(ACT/'mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits')),nside)**2
    elif survey=='Planck':
        # the PR4 products are in GALACTIC coordinates; the tracers and ACT are equatorial (celestial)
        alm=np.nan_to_num(hp.read_alm(str(PLANCK/'PR42018like_klm_dat_MV.fits')))
        rot=hp.Rotator(coord=['G','C']); alm=rot.rotate_alm(alm)
        mask=rot.rotate_map_pixel(hp.ud_grade(hp.read_map(str(PLANCK/'mask.fits.gz')),nside))
        mask=(mask>0.5).astype(float)      # binary again after the pixel rotation
    else: raise ValueError(survey)
    lm=hp.Alm.getlmax(len(alm)); ell=np.arange(lm+1)
    alm=hp.almxfl(alm,((ell>=2)&(ell<=lmax)).astype(float))
    return hp.alm2map(alm,nside,verbose=False),np.asarray(mask,float)
