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
import numpy as np
import h5py
import healpy as hp
from astropy.io import fits
from lyalenser.paths import DESI_DELTAS
from lyalenser.cosmo import chi as chi_of_z
from lyalenser.config import Config, SightlineSet

LYA=1215.67
LYB=1025.72
# The VAC ships the Lya forest in two rest-frame windows: region A (1040-1205 A, the Lya forest proper) and
# region B (920-1020 A, blueward of Lyb emission).  A region-B pixel carries Lya absorption at
# z = lambda/1215.67 - 1, which is the field we want, AND Lyb absorption at z = lambda/1025.72 - 1, which is a
# different and much more distant slab.  For pixel pairs closer than 30 Mpc/h the Lya-Lyb cross term connects
# two fields separated by ~1000 Mpc/h radially and is negligible, so A x B pairs measure Lya-Lya alone; B x B
# pairs also carry Lyb-Lyb at the same z and are dropped (see `pairs.accumulate`).
REGION_DIRS={'lya':'delta-lya-0-0','lyb':'delta-lyb-0-0'}
REGION_CODE={'lya':0,'lyb':1}


def delta_files(region='lya'):
    d=DESI_DELTAS if region=='lya' else DESI_DELTAS.parent.parent/REGION_DIRS[region]/'Delta'
    return sorted(glob.glob(str(d/'delta-*.fits.gz')))


def in_region(ra_deg,dec_deg,region):
    if region is None: return np.ones(len(ra_deg),bool)
    if isinstance(region,dict) and 'disc' in region:
        ra0,dec0,rad=region['disc']
        v=hp.ang2vec(ra_deg,dec_deg,lonlat=True); v0=hp.ang2vec(ra0,dec0,lonlat=True)
        return np.degrees(np.arccos(np.clip(v@v0,-1,1)))<=rad
    if isinstance(region,dict) and 'nside' in region:
        pix=hp.ang2pix(region['nside'],ra_deg,dec_deg,lonlat=True); return np.isin(pix,region['pixels'])
    raise ValueError('region must be {"disc": (ra, dec, radius_deg)} or {"nside": n, "pixels": [...]}')


def _scan_region(region,zmin,zmax,sky,min_pixels,slab_edges,files=None,max_files=None,verbose=True):
    """Per-forest records of one delta region, selected on the Lya absorption redshift."""
    files=delta_files(region) if files is None else list(files)
    if max_files: files=files[:max_files]
    out=[]; t0=time.perf_counter(); nread=0
    for i,f in enumerate(files):
        with fits.open(f,memmap=False) as h:
            md=h['METADATA'].data; rad=np.degrees(md['RA']); decd=np.degrees(md['DEC'])
            keep=in_region(rad,decd,sky)
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
                out.append(dict(qid=int(md['LOS_ID'][keep][j]),ra=float(rad[keep][j]),dec=float(decd[keep][j]),
                                zq=float(md['Z'][keep][j]),snr=float(md['MEANSNR'][keep][j]),chi=c,
                                delta=d[j][m].astype(np.float32),w=w[j][m].astype(np.float32),slab=lab,
                                region=np.full(len(c),REGION_CODE[region],np.int8)))
            nread+=1
        if verbose and i%100==0: print(f'  [{region}] {i}/{len(files)} files, {len(out)} forests, {time.perf_counter()-t0:.0f} s',flush=True)
    if verbose: print(f'  [{region}] {len(out)} forests from {nread} files, {time.perf_counter()-t0:.0f} s',flush=True)
    return out


def read_deltas(zmin=2.1,zmax=3.0,region=None,min_pixels=50,cfg=None,files=None,max_files=None,verbose=True,
                forest_regions=('lya',)):
    """SightlineSet over the requested delta regions.

    ``forest_regions`` selects the rest-frame windows: ('lya',) reproduces the region-A measurement exactly,
    ('lya','lyb') adds the Lyb-region segment of every quasar that has one as an EXTENSION OF THE SAME SIGHTLINE
    (the pixels are concatenated and sorted by distance, and each carries a region tag).  The >= ``min_pixels``
    cut is applied per region, so adding region B leaves the region-A sample bit-identical.
    """
    cfg=Config() if cfg is None else cfg
    slab_edges=[(float(chi_of_z(a)),float(chi_of_z(b))) for a,b in cfg.slabs]
    t0=time.perf_counter()
    recs=[]
    for r in forest_regions:
        recs+= _scan_region(r,zmin,zmax,region,min_pixels,slab_edges,
                            files=files if r=='lya' else None,max_files=max_files,verbose=verbose)
    if not recs: raise RuntimeError('no forests selected')
    by_q={}
    for rec in recs: by_q.setdefault(rec['qid'],[]).append(rec)
    qid=[]; ra=[]; dec=[]; zq=[]; chis=[]; deltas=[]; weights=[]; slabs=[]; regions=[]; npix=[]; snr=[]
    for q,group in by_q.items():
        group.sort(key=lambda g: g['region'][0])            # region A first, so its metadata wins
        c=np.concatenate([g['chi'] for g in group])
        order=np.argsort(c,kind='stable')
        qid.append(q); ra.append(group[0]['ra']); dec.append(group[0]['dec']); zq.append(group[0]['zq'])
        snr.append(group[0]['snr'])
        chis.append(c[order])
        for key,acc in (('delta',deltas),('w',weights),('slab',slabs),('region',regions)):
            acc.append(np.concatenate([g[key] for g in group])[order])
        npix.append(len(order))
    starts=np.r_[0,np.cumsum(npix)].astype(np.int64)
    attrs={'zmin':float(zmin),'zmax':float(zmax),'description':'DESI DR1 Lya deltas (picca, DELTA_BLIND, WEIGHT)','chi_ref':float(cfg.chi_ref),
           'source':'+'.join(REGION_DIRS[r] for r in forest_regions),'min_pixels':int(min_pixels),'region':str(region),
           'forest_regions':'+'.join(forest_regions),
           'weights':'picca WEIGHT (inverse variance, LSS variance model; the C^-1 diagonal)','pixel_A':0.8}
    sl=SightlineSet(np.asarray(qid),np.asarray(ra),np.asarray(dec),np.asarray(zq,np.float32),starts,np.concatenate(chis),
                    np.concatenate(deltas),np.concatenate(weights),np.concatenate(slabs),attrs,
                    np.concatenate(regions))
    sl.attrs['meansnr']=np.asarray(snr,np.float32)
    if verbose:
        nb=int((sl.region==1).sum())
        print(f'read {sl.nq} sightlines, {len(sl.chi)} pixels ({nb} in region B, {100*nb/len(sl.chi):.1f} %), '
              f'{time.perf_counter()-t0:.0f} s',flush=True)
    return sl


def save_sightlines(sl,path):
    import h5py
    with h5py.File(path,'w') as f:
        g=f.create_group('sightlines')
        for k in ('qid','ra','dec','zq','pix_start','chi','delta','w','slab','region'): g[k]=getattr(sl,k)
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


def load_sightlines(path):
    """Read a `SightlineSet` written by `save_sightlines` (group `sightlines` of an HDF5 file)."""
    with h5py.File(path, "r") as f:
        g = f["sightlines"]; vals = [g[k][()] for k in ("qid", "ra", "dec", "zq", "pix_start", "chi", "delta", "w", "slab")]
        region = g["region"][()] if "region" in g else None
        attrs = dict(g.attrs)
    return SightlineSet(*vals, attrs, region)
