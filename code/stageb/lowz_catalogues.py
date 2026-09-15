"""DESI DR1 low-redshift tracers -> per-slice kappa_lya templates on the sphere (Stage B, iteration 7).

For every (tracer, redshift slice) of lowz.TRACERS: read the LSS clustering catalogue (NGC + SGC) and its randoms,
select the slice, build the kernel-weighted map (forest-source lensing kernel W(chi; chi_ref), per-object weight
WEIGHT / (nbar(chi) Omega_pix), randoms subtracted, completeness from the randoms), fit the linear bias from the
masked pseudo-C_ell on large scales (40 <= ell <= 0.2 chi(z_mid); shot noise subtracted; theory C_ll(slice) x
pixel window, fsky-scaled), and Wiener-combine the tracers of a slice with the model covariance (theory signal
+ shot noise) into the slice's kappa_lya estimate; the combined estimate is the sum over slices. Outputs go to
$LYALENSER_DATA/lowz/: unit-bias maps, masks, filtered alm per slice and combined, and a JSON summary.
The pair estimator on DR1 sightlines (Stage B proper) consumes the alm through templates.alpha_at.

Usage: python lowz_catalogues.py [--nside 512] [--tracers LRG ELG QSO] [--out DIR]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
import fitsio
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline'):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import RAW, DATA
from cosmo import chi as chi_of_z, z_of_chi
from cross_spectrum import kernel
from config import Config
from templates import matched_template
from lowz import TRACERS, slices_of, slice_spectra, bias_band, fit_bias, ANNULUS, Tracer

LSS=RAW/'desi/lss_v1.5'
FILES={'LRG':'LRG','ELG':'ELG_LOPnotqso','QSO':'QSO'}
RANDOMS_PER_CAP=2


def read_catalogue(tracer,zmin,zmax,randoms=RANDOMS_PER_CAP):
    """Data and randoms of both caps in the slice: ra, dec, z, weight."""
    def read(path):
        d=fitsio.read(str(path),columns=['RA','DEC','Z','WEIGHT'])
        keep=(d['Z']>=zmin)&(d['Z']<zmax)
        return {k.lower():np.asarray(d[k][keep],float) for k in ('RA','DEC','Z','WEIGHT')}
    cat={k:[] for k in ('ra','dec','z','weight')}; rnd={k:[] for k in ('ra','dec','z','weight')}
    for cap in ('NGC','SGC'):
        d=read(LSS/f'{FILES[tracer]}_{cap}_clustering.dat.fits')
        for k in cat: cat[k].append(d[k])
        for i in range(randoms):
            r=read(LSS/f'{FILES[tracer]}_{cap}_{i}_clustering.ran.fits')
            for k in rnd: rnd[k].append(r[k])
    return {k:np.concatenate(v) for k,v in cat.items()},{k:np.concatenate(v) for k,v in rnd.items()}


def footprint(rnd,nside,frac=.5):
    """Pixels holding at least `frac` of the mean random count over the survey (the mean is taken over pixels with
    any randoms): partially covered rim pixels are excluded, otherwise they dilute the fsky-scaled pseudo-C_ell
    (a 1-random threshold gave fsky 0.185 for the LRGs against ~0.14 of actual footprint, i.e. 25 % low power)."""
    pix=hp.ang2pix(nside,rnd['ra'],rnd['dec'],lonlat=True); c=np.bincount(pix,minlength=hp.nside2npix(nside))
    return c>=frac*c[c>0].mean()


def binned_cl(cl,lmax,edges):
    ell=np.arange(lmax+1); L=.5*(edges[1:]+edges[:-1]); out=np.zeros(len(L)); nm=np.zeros(len(L))
    for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        m=(ell>=a)&(ell<b); w=2*ell[m]+1
        if w.sum()>0: out[i]=np.sum(w*cl[m])/w.sum(); nm[i]=w.sum()
    return L,out,nm


def mask_coupled_sphere(theory_cl,mask,lmax,edges,nsim=16,seed=0):
    """Expected binned pseudo-C_ell of a Gaussian field with spectrum theory_cl (already including the pixel window)
    observed through `mask`, divided by fsky = <mask^2>: Monte-Carlo transfer instead of the fsky approximation
    (the DR1 footprint is fragmented, and the fsky scaling biases low-ell bandpowers)."""
    nside=hp.npix2nside(len(mask)); rng=np.random.default_rng(seed); fsky=float(np.mean(mask**2)); acc=None
    state=np.random.get_state()
    try:
        for i in range(nsim):
            np.random.seed(int(rng.integers(2**31)))
            m=hp.synfast(theory_cl,nside,lmax=lmax,pixwin=False,verbose=False)
            cl=hp.anafast(m*mask,lmax=lmax,iter=0)/max(fsky,1e-30)
            _,b,_=binned_cl(cl,lmax,edges); acc=b if acc is None else acc+b
    finally: np.random.set_state(state)
    return acc/nsim


def build_slice(slice_,tracers,cfg,nside,lmax,out):
    """Unit-bias maps, bias fits and the Wiener-combined kappa_lya estimate of one slice."""
    cref=cfg.chi_ref; edges=np.arange(0,lmax+ANNULUS,ANNULUS); L=.5*(edges[1:]+edges[:-1])
    Lth,C=slice_spectra(slice_['zmin'],slice_['zmax'],cref,lmax); S=np.interp(np.arange(lmax+1),Lth,C[:,1,1])
    pw=hp.pixwin(nside,lmax=lmax); maps=[]; shots=[]; info={'zmin':slice_['zmin'],'zmax':slice_['zmax'],'tracers':{}}
    alms=[]; masks=[]
    for t in tracers:
        t0=time.perf_counter(); cat,rnd=read_catalogue(t.name,t.zmin,t.zmax); fp=footprint(rnd,nside)
        rpix=hp.ang2pix(nside,rnd['ra'],rnd['dec'],lonlat=True); rc=np.bincount(rpix,minlength=hp.nside2npix(nside)); nmin=int(.5*rc[rc>0].mean())
        unit=lambda z: np.ones_like(np.asarray(z,float))
        def template(c,r):
            # the random-count threshold scales with the random set actually used (halves have half the randoms)
            rp=hp.ang2pix(nside,r['ra'],r['dec'],lonlat=True); rcnt=np.bincount(rp,minlength=hp.nside2npix(nside))
            return matched_template({'ra':c['ra'],'dec':c['dec'],'z':c['z']},{'ra':r['ra'],'dec':r['dec'],'z':r['z']},unit,cfg,footprint_mask=fp,
                                    source_chi=cref,radial_bins=10,data_weights=c['weight'],random_weights=r['weight'],nside=nside,lmax=lmax,
                                    nmin_rand=int(.5*rcnt[rcnt>0].mean()))
        alm,mask,meta=template(cat,rnd)
        # Split estimator: two random halves of the data (and of the randoms) give a cross-spectrum free of shot
        # noise (the bias) and a half-difference whose power is the shot noise of the full map (the Wiener model),
        # including every weight, completeness and mask effect without a model.
        rng=np.random.default_rng(12345); hd=rng.random(len(cat['ra']))<.5; hr=rng.random(len(rnd['ra']))<.5
        half=lambda d,m: {k:v[m] for k,v in d.items()}
        _,maskA,ma=template(half(cat,hd),half(rnd,hr)); _,maskB,mb=template(half(cat,~hd),half(rnd,~hr))
        mask=mask&maskA&maskB   # common mask of the full map and both halves
        mA=ma['kappa_map']*mask; mB=mb['kappa_map']*mask; fsky=float(np.mean(mask.astype(float)**2)); kmap=meta['kappa_map']*mask
        alm=hp.map2alm(kmap,lmax=lmax,iter=0)
        cross=hp.anafast(mA,mB,lmax=lmax,iter=0)/max(fsky,1e-30); autoA=hp.anafast(mA,lmax=lmax,iter=0)/fsky; autoB=hp.anafast(mB,lmax=lmax,iter=0)/fsky
        diff=hp.anafast(.5*(mA-mB),lmax=lmax,iter=0)/fsky; full=hp.anafast(kmap*mask,lmax=lmax,iter=0)/fsky
        Lb,meas,nm=binned_cl(cross,lmax,edges); _,aA,_=binned_cl(autoA,lmax,edges); _,aB,_=binned_cl(autoB,lmax,edges); _,dd,_=binned_cl(diff,lmax,edges); _,fa,_=binned_cl(full,lmax,edges)
        _,Tfsky,_=binned_cl(S*pw**2,lmax,edges); T=mask_coupled_sphere(S*pw**2,mask.astype(float),lmax,edges)
        band=(Lb>=bias_band(t.zmin,t.zmax)[0])&(Lb<=bias_band(t.zmin,t.zmax)[1])&(T>0)&(nm>0)
        var=(aA*aB+meas**2)/np.maximum(nm*fsky,1)          # Gaussian variance of the cross bandpower
        w=1/np.maximum(var,1e-40); b2=float(np.sum((w*meas*T)[band])/np.sum((w*T*T)[band])); vb2=1/float(np.sum((w*T*T)[band]))
        b=float(np.sqrt(max(b2,1e-12))); sb=float(np.sqrt(vb2)/(2*b)); chi2=float(np.sum((w*(meas-b2*T)**2)[band]))
        shot=float(np.mean(dd[(Lb>=100)&(Lb<=lmax-40)]))    # white level of the half-difference = shot noise of the full map
        fit={'b':b,'sigma_b':sb,'b2':b2,'sigma_b2':float(np.sqrt(vb2)),'chi2':chi2,'dof':int(band.sum()-1),'band':[float(x) for x in bias_band(t.zmin,t.zmax)],
             'annuli_used':int(band.sum()),'estimator':'cross-spectrum of two random halves (no shot noise)',
             'mask_transfer_band':(T/np.maximum(Tfsky,1e-30))[(Lb>=40)&(Lb<=300)].tolist(),'shot_from_half_difference':shot,'shot_model':meta['shot_s']}
        usable=np.isfinite(fit['b2']) and fit['b2']>3*fit['sigma_b2']; b=fit['b'] if usable else t.bias
        fit['usable']=bool(usable); fit['bias_source']='auto-spectrum' if usable else 'fallback_table'
        maps.append(kmap/b); alms.append(alm/b); masks.append(mask); shots.append(shot/b**2)
        info['tracers'][t.label]={'n_objects':int(len(cat['ra'])),'sum_weights':float(cat['weight'].sum()),'fsky':fsky,'bias_fit':fit,'bias_used':float(b),
                                  'bias_table':t.bias,'shot_s_unit_bias':shot,'shot_s_uniform_model':meta['shot_s_uniform_model'],
                                  'cross_cl_binned':meas.tolist(),'full_auto_binned':fa.tolist(),'half_difference_binned':dd.tolist(),
                                  'theory_unit_bias_binned':T.tolist(),'L':Lb.tolist(),'wall_s':time.perf_counter()-t0}
        hp.write_map(str(out/f'unitbias_{t.label}_nside{nside}.fits'),[kmap,mask.astype(float)],overwrite=True,dtype=np.float64)
        print(f"  {t.label}: {len(cat['ra'])} objects, fsky {fsky:.3f}, b = {fit['b']:.3f} +- {fit['sigma_b']:.3f} ({fit['bias_source']}; shot half-diff/model {shot/meta['shot_s']:.2f}), {time.perf_counter()-t0:.0f} s",flush=True)
    # Wiener combination with the model covariance C_kl = S pw^2 + shot_k delta_kl, signal S pw (cross with the continuous field).
    k=len(maps); mask=np.prod(masks,axis=0).astype(bool); ell=np.arange(lmax+1); W=np.zeros((lmax+1,k))
    for l in range(2,lmax+1):
        Cm=S[l]*pw[l]**2*np.ones((k,k))+np.diag(shots); W[l]=np.linalg.solve(Cm,np.full(k,S[l]*pw[l]))
    comb=sum(hp.almxfl(alms[i],W[:,i]) for i in range(k))
    hp.write_alm(str(out/f"kappa_slice_{slice_['zmin']:g}_{slice_['zmax']:g}_alm.fits"),comb,overwrite=True)
    hp.write_map(str(out/f"mask_slice_{slice_['zmin']:g}_{slice_['zmax']:g}_nside{nside}.fits"),mask.astype(float),overwrite=True,dtype=np.float64)
    info['weights_at_L']={str(l):W[l].tolist() for l in (40,100,200,300)}; info['shot_s']=shots; info['fsky']=float(mask.mean())
    return comb,mask,info


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--nside',type=int,default=512); ap.add_argument('--lmax',type=int,default=1000)
    ap.add_argument('--tracers',nargs='*',default=['LRG','ELG','QSO']); ap.add_argument('--out',type=Path,default=DATA/'lowz')
    ap.add_argument('--slices',nargs='*',type=float,default=None,help='zmin zmax pairs to restrict to')
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True); cfg=Config(scale=1.,r_perp_min=3.,fit_rperp_min=3.)
    tracers=[t for t in TRACERS if t.name in a.tracers]; summary={'nside':a.nside,'lmax':a.lmax,'slices':[]}; total=None; total_mask=None
    for s in slices_of(tuple(tracers)):
        if a.slices and not any(abs(s['zmin']-a.slices[i])<1e-6 and abs(s['zmax']-a.slices[i+1])<1e-6 for i in range(0,len(a.slices),2)): continue
        print(f"slice {s['zmin']}-{s['zmax']}: {[t.label for t in s['tracers']]}",flush=True)
        comb,mask,info=build_slice(s,s['tracers'],cfg,a.nside,a.lmax,a.out); summary['slices'].append(info)
        total=comb if total is None else total+comb; total_mask=mask if total_mask is None else (total_mask&mask)
        (a.out/'summary.json').write_text(json.dumps(summary,indent=1,default=float)+'\n')
    if total is not None:
        hp.write_alm(str(a.out/'kappa_combined_alm.fits'),total,overwrite=True)
        hp.write_map(str(a.out/f'mask_combined_nside{a.nside}.fits'),total_mask.astype(float),overwrite=True,dtype=np.float64)
    print('done',a.out)


if __name__=='__main__': main()
