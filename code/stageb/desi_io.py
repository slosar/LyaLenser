"""DESI DR1 Lya deltas -> SightlineSet on the sphere (Stage B).

Delta files (picca, `delta-lya-0-0/Delta/delta-*.fits.gz`): HDU LAMBDA (2716 wavelengths, 0.8 A linear grid,
3600-5772 A), METADATA (LOS_ID, RA, DEC in radians, Z, MEANSNR, TARGETID), DELTA_BLIND (forests x wavelengths;
the blinding only affects picca's fiducial distance conversion, the arrays are the measured fluctuations),
WEIGHT (picca inverse-variance weights with the LSS variance model: the C^-1 diagonal of the pair estimator),
CONT. Masked or absent pixels have WEIGHT = 0. DLAs and BALs are masked upstream by the DR1 delta extraction.

`read_deltas(zmin, zmax, region=...)` keeps pixels with zmin <= z <= zmax, forests with >= min_pixels kept
pixels inside an optional sky region (HEALPix pixel list or an (ra, dec, radius) disc), and returns the
SightlineSet (ra, dec in degrees; chi in Mpc/h; per-pixel weights; slab labels from cfg.slabs).
"""
from __future__ import annotations
import glob, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
from astropy.io import fits
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline'):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DESI_DELTAS
from cosmo import chi as chi_of_z
from config import Config, SightlineSet

LYA=1215.67


def delta_files():
    return sorted(glob.glob(str(DESI_DELTAS/'delta-*.fits.gz')))


def in_region(ra_deg,dec_deg,region):
    if region is None: return np.ones(len(ra_deg),bool)
    if isinstance(region,dict) and 'disc' in region:
        ra0,dec0,rad=region['disc']
        v=hp.ang2vec(ra_deg,dec_deg,lonlat=True); v0=hp.ang2vec(ra0,dec0,lonlat=True)
        return np.degrees(np.arccos(np.clip(v@v0,-1,1)))<=rad
    if isinstance(region,dict) and 'nside' in region:
        pix=hp.ang2pix(region['nside'],ra_deg,dec_deg,lonlat=True); return np.isin(pix,region['pixels'])
    raise ValueError('region must be {"disc": (ra, dec, radius_deg)} or {"nside": n, "pixels": [...]}')


def read_deltas(zmin=2.1,zmax=3.0,region=None,min_pixels=50,cfg=None,files=None,max_files=None,verbose=True):
    cfg=Config() if cfg is None else cfg
    files=delta_files() if files is None else list(files)
    if max_files: files=files[:max_files]
    slab_edges=[(float(chi_of_z(a)),float(chi_of_z(b))) for a,b in cfg.slabs]
    qid=[]; ra=[]; dec=[]; zq=[]; chis=[]; deltas=[]; weights=[]; slabs=[]; npix=[]; snr=[]
    t0=time.perf_counter(); nread=0
    for i,f in enumerate(files):
        with fits.open(f,memmap=False) as h:
            md=h['METADATA'].data; rad=np.degrees(md['RA']); decd=np.degrees(md['DEC'])
            keep=in_region(rad,decd,region)
            if not keep.any(): continue
            lam=h['LAMBDA'].data; z=lam/LYA-1; inz=(z>=zmin)&(z<=zmax)
            if not inz.any(): continue
            chi_l=np.asarray(chi_of_z(z[inz]),np.float32)
            d=h['DELTA_BLIND'].data[keep][:,inz]; w=h['WEIGHT'].data[keep][:,inz]
            ok=np.isfinite(d)&(w>0)
            for j in np.flatnonzero(ok.sum(axis=1)>=min_pixels):
                m=ok[j]; c=chi_l[m]
                lab=np.full(len(c),-1,np.int8)
                for k,(lo,hi) in enumerate(slab_edges): lab[(c>=lo)&(c<hi)]=k
                qid.append(int(md['LOS_ID'][keep][j])); ra.append(float(rad[keep][j])); dec.append(float(decd[keep][j])); zq.append(float(md['Z'][keep][j]))
                snr.append(float(md['MEANSNR'][keep][j]))
                chis.append(c); deltas.append(d[j][m].astype(np.float32)); weights.append(w[j][m].astype(np.float32)); slabs.append(lab); npix.append(int(m.sum()))
            nread+=1
        if verbose and i%100==0: print(f'  {i}/{len(files)} files, {len(qid)} forests, {time.perf_counter()-t0:.0f} s',flush=True)
    if not qid: raise RuntimeError('no forests selected')
    starts=np.r_[0,np.cumsum(npix)].astype(np.int64)
    attrs={'zmin':float(zmin),'zmax':float(zmax),'description':'DESI DR1 Lya deltas (picca, DELTA_BLIND, WEIGHT)','chi_ref':float(cfg.chi_ref),
           'source':'delta-lya-0-0','min_pixels':int(min_pixels),'n_files_used':int(nread),'region':str(region),
           'weights':'picca WEIGHT (inverse variance, LSS variance model; the C^-1 diagonal)','pixel_A':0.8}
    sl=SightlineSet(np.asarray(qid),np.asarray(ra),np.asarray(dec),np.asarray(zq,np.float32),starts,np.concatenate(chis),
                    np.concatenate(deltas),np.concatenate(weights),np.concatenate(slabs),attrs)
    sl.attrs['meansnr']=np.asarray(snr,np.float32)
    if verbose: print(f'read {sl.nq} forests, {len(sl.chi)} pixels, {time.perf_counter()-t0:.0f} s',flush=True)
    return sl


def save_sightlines(sl,path):
    import h5py
    with h5py.File(path,'w') as f:
        g=f.create_group('sightlines')
        for k in ('qid','ra','dec','zq','pix_start','chi','delta','w','slab'): g[k]=getattr(sl,k)
        for k,v in sl.attrs.items():
            if isinstance(v,np.ndarray): g[f'attr_{k}']=v
            else: g.attrs[k]=v


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--disc',type=float,nargs=3,metavar=('RA','DEC','RADIUS_DEG'))
    ap.add_argument('--zmin',type=float,default=2.1); ap.add_argument('--zmax',type=float,default=3.0); ap.add_argument('--out')
    a=ap.parse_args()
    sl=read_deltas(a.zmin,a.zmax,region={'disc':tuple(a.disc)} if a.disc else None)
    print('forests',sl.nq,'pixels',len(sl.chi),'median pixels per forest',np.median(np.diff(sl.pix_start)),'chi range',sl.chi.min(),sl.chi.max())
    if a.out: save_sightlines(sl,a.out)
