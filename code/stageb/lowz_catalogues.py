"""Low-redshift tracers -> per-slice kappa_lya templates on the sphere (Stage B, iteration 10).

For every (tracer, redshift slice) of TRACERS: read the clustering catalogue and its randoms (DESI DR1 LSS v1.5:
LRG, ELG_LOPnotqso, QSO, BGS_BRIGHT-21.5; BOSS DR12v5 CMASSLOWZTOT), select the slice, build the kernel-weighted
map (forest-source lensing kernel W(chi; chi_ref), per-object weight WEIGHT / (nbar(chi) Omega_pix), randoms
subtracted, completeness from the randoms), and fit the linear bias from the NaMaster cross-spectrum of two random
halves of the catalogue (no shot noise) against the Limber spectrum of the slice pushed through the same bandpower
windows, on large scales only (40 <= ell <= 0.2 chi(z_mid)), with the NaMaster Gaussian covariance.

Tracers of a slice are Wiener-combined into the slice's kappa_lya estimate with the model covariance
C_kl = S pw^2 + N_kl (S the slice's convergence spectrum, N the shot noise: diagonal from the half-difference
spectrum of each map, off-diagonal from the objects two catalogues share, e.g. BOSS CMASS galaxies that are also
DESI LRGs). Because the tracer footprints differ (BOSS covers sky DR1 does not, and vice versa), the weights are
computed PER COVERAGE CLASS (the subset of tracers covering a pixel): within each class the filtered sum is the
conditional expectation of kappa_lya given the tracers actually present, which keeps the amplitude estimator
normalised to A = 1 without further calibration. Slices are independent, so the combined estimate is the sum over
slices; the combined mask is the union of the slice masks. Outputs go to --out: unit-bias maps + masks, filtered
alm per slice and combined, and summary.json (biases, spectra, shot noise, class weights).

Usage: python lowz_catalogues.py [--nside 512] [--lmax 1000] [--tracers LRG ELG QSO BGS BOSS] [--out DIR]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
from itertools import combinations
import numpy as np
import healpy as hp
import fitsio
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import RAW, DATA
from cosmo import chi as chi_of_z, z_of_chi
from cross_spectrum import kernel
from config import Config
from templates import matched_template
from lowz import slices_of, slice_spectra, bias_band, ANNULUS, Tracer
from nmt_spectra import Spectra, fit_amplitude

LSS=RAW/'desi/lss_v1.5'
BOSS=RAW/'boss'
DESI_FILES={'LRG':'LRG','ELG':'ELG_LOPnotqso','QSO':'QSO','BGS':'BGS_BRIGHT-21.5'}
RANDOMS_PER_CAP=2
# (name, zmin, zmax, fallback bias). The bias is measured from the auto-spectrum; the table value is only used when
# the measurement is not significant. BOSS enters three slices (LOWZ below 0.43, CMASS above), BGS the lowest one.
TRACERS=(Tracer('BGS',.1,.4,1.4,0.),Tracer('BOSS',.1,.4,1.8,0.),
         Tracer('LRG',.4,.6,1.9,0.),Tracer('BOSS',.4,.6,2.0,0.),
         Tracer('LRG',.6,.8,2.1,0.),Tracer('BOSS',.6,.8,2.1,0.),
         Tracer('LRG',.8,1.1,2.3,0.),Tracer('ELG',.8,1.1,1.3,0.),Tracer('QSO',.8,1.1,1.7,0.),
         Tracer('ELG',1.1,1.6,1.4,0.),Tracer('QSO',1.1,1.6,2.1,0.),
         Tracer('QSO',1.6,1.75,2.5,0.))
MATCH_ARCSEC=1.0      # positional match of shared objects between two catalogues of the same slice


def _select(d,zmin,zmax):
    keep=(d['z']>=zmin)&(d['z']<zmax)
    return {k:v[keep] for k,v in d.items()}


def read_catalogue(tracer,zmin,zmax,randoms=RANDOMS_PER_CAP):
    """Data and randoms of both caps in the slice: ra, dec, z, weight (DESI: WEIGHT; BOSS: SYSTOT (CP + NOZ - 1),
    randoms unweighted)."""
    def desi(path):
        d=fitsio.read(str(path),columns=['RA','DEC','Z','WEIGHT'])
        return _select({k.lower():np.asarray(d[k],float) for k in ('RA','DEC','Z','WEIGHT')},zmin,zmax)
    def boss(path,is_random):
        cols=['RA','DEC','Z'] if is_random else ['RA','DEC','Z','WEIGHT_SYSTOT','WEIGHT_CP','WEIGHT_NOZ']
        d=fitsio.read(str(path),columns=cols)
        out={k.lower():np.asarray(d[k],float) for k in ('RA','DEC','Z')}
        out['weight']=np.ones(len(out['ra'])) if is_random else np.asarray(d['WEIGHT_SYSTOT']*(d['WEIGHT_CP']+d['WEIGHT_NOZ']-1.),float)
        return _select(out,zmin,zmax)
    cat={k:[] for k in ('ra','dec','z','weight')}; rnd={k:[] for k in ('ra','dec','z','weight')}
    if tracer=='BOSS':
        for cap in ('North','South'):
            d=boss(BOSS/f'galaxy_DR12v5_CMASSLOWZTOT_{cap}.fits.gz',False)
            for k in cat: cat[k].append(d[k])
            r=boss(BOSS/f'random0_DR12v5_CMASSLOWZTOT_{cap}.fits.gz',True)
            for k in rnd: rnd[k].append(r[k])
    else:
        for cap in ('NGC','SGC'):
            d=desi(LSS/f'{DESI_FILES[tracer]}_{cap}_clustering.dat.fits')
            for k in cat: cat[k].append(d[k])
            for i in range(randoms):
                r=desi(LSS/f'{DESI_FILES[tracer]}_{cap}_{i}_clustering.ran.fits')
                for k in rnd: rnd[k].append(r[k])
    return {k:np.concatenate(v) for k,v in cat.items()},{k:np.concatenate(v) for k,v in rnd.items()}


def footprint(rnd,nside,frac=.5):
    """Pixels holding at least `frac` of the mean random count over the survey (the mean is taken over pixels with
    any randoms): partially covered rim pixels are excluded."""
    pix=hp.ang2pix(nside,rnd['ra'],rnd['dec'],lonlat=True); c=np.bincount(pix,minlength=hp.nside2npix(nside))
    return c>=frac*c[c>0].mean()


def binned_cl(cl,lmax,edges):
    ell=np.arange(lmax+1); L=.5*(edges[1:]+edges[:-1]); out=np.zeros(len(L)); nm=np.zeros(len(L))
    for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        m=(ell>=a)&(ell<b); w=2*ell[m]+1
        if w.sum()>0: out[i]=np.sum(w*cl[m])/w.sum(); nm[i]=w.sum()
    return L,out,nm


def shared_objects(c1,c2,arcsec=MATCH_ARCSEC,dz=0.01):
    """Index pairs (i1, i2) of objects present in both catalogues (position within ``arcsec``, redshift within
    ``dz`` (1 + z))."""
    from scipy.spatial import cKDTree
    def vec(c): r=np.deg2rad(c['ra']); d=np.deg2rad(c['dec']); return np.column_stack((np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)))
    if len(c1['ra'])==0 or len(c2['ra'])==0: return np.zeros(0,int),np.zeros(0,int)
    t=cKDTree(vec(c2)); dist,j=t.query(vec(c1),distance_upper_bound=2*np.sin(np.deg2rad(arcsec/3600)/2))
    ok=np.isfinite(dist)
    i=np.flatnonzero(ok); j=j[ok]
    same=np.abs(c1['z'][i]-c2['z'][j])<dz*(1+c1['z'][i])
    return i[same],j[same]


def _unit_weights(cat,rnd,fp,cfg,nside,cref):
    """Per-object map weights u_i w_i of `matched_template` (needed for the cross shot noise of shared objects)."""
    npix=hp.nside2npix(nside); area=4*np.pi/npix; omega=np.count_nonzero(fp)*area
    qc=chi_of_z(cat['z']); c1=min(float(qc.min()),float(chi_of_z(rnd['z']).min())); c2=max(float(qc.max()),float(chi_of_z(rnd['z']).max()))
    edges=np.linspace(c1,c2,11); hist,_=np.histogram(qc,edges,weights=cat['weight'])
    ib=np.clip(np.searchsorted(edges,qc,side='right')-1,0,len(hist)-1)
    nbar=np.maximum(hist[ib]/(np.diff(edges)[ib]*max(omega,area)),1e-30)
    return kernel(qc,cref)/(nbar*area)*cat['weight']


def coverage_classes(masks):
    """Non-empty subsets of tracers as (tuple of indices, class mask)."""
    k=len(masks); out=[]
    for r in range(1,k+1):
        for sub in combinations(range(k),r):
            m=np.ones(len(masks[0]),bool)
            for i in range(k): m&=masks[i] if i in sub else ~masks[i]
            if m.any(): out.append((sub,m))
    return out


def build_slice(slice_,tracers,cfg,nside,lmax,out,S:Spectra):
    """Unit-bias maps, bias fits and the Wiener-combined kappa_lya estimate of one slice."""
    cref=cfg.chi_ref
    Lth,C=slice_spectra(slice_['zmin'],slice_['zmax'],cref,lmax); Sl=np.interp(np.arange(lmax+1),Lth,C[:,1,1])
    pw=hp.pixwin(nside,lmax=lmax); ell=np.arange(lmax+1)
    info={'zmin':slice_['zmin'],'zmax':slice_['zmax'],'tracers':{},'estimator':'NaMaster (pymaster) decoupled bandpowers, width 40'}
    maps=[]; masks=[]; shots=[]; cats=[]; uw=[]; biases=[]
    for t in tracers:
        t0=time.perf_counter(); cat,rnd=read_catalogue(t.name,t.zmin,t.zmax); fp=footprint(rnd,nside)
        unit=lambda z: np.ones_like(np.asarray(z,float))
        def template(c,r):
            rp=hp.ang2pix(nside,r['ra'],r['dec'],lonlat=True); rcnt=np.bincount(rp,minlength=hp.nside2npix(nside))
            return matched_template({'ra':c['ra'],'dec':c['dec'],'z':c['z']},{'ra':r['ra'],'dec':r['dec'],'z':r['z']},unit,cfg,footprint_mask=fp,
                                    source_chi=cref,radial_bins=10,data_weights=c['weight'],random_weights=r['weight'],nside=nside,lmax=lmax,
                                    nmin_rand=int(.5*rcnt[rcnt>0].mean()))
        _,mask,meta=template(cat,rnd)
        # Split estimator: two random halves of the data (and of the randoms) give a cross-spectrum free of shot
        # noise (the bias) and a half-difference whose power is the shot noise of the full map (the Wiener model).
        rng=np.random.default_rng(12345); hd=rng.random(len(cat['ra']))<.5; hr=rng.random(len(rnd['ra']))<.5
        half=lambda d,m: {k:v[m] for k,v in d.items()}
        _,maskA,ma=template(half(cat,hd),half(rnd,hr)); _,maskB,mb=template(half(cat,~hd),half(rnd,~hr))
        mask=mask&maskA&maskB; fsky=float(mask.mean())
        mA=ma['kappa_map']*mask; mB=mb['kappa_map']*mask; kmap=meta['kappa_map']*mask; diff=.5*(mA-mB)
        key=f'{t.label}_mask'
        fA=S.field(mask,[mA],key=key); fB=S.field(mask,[mB],key=key); fF=S.field(mask,[kmap],key=key); fD=S.field(mask,[diff],key=key)
        cross=S.cross(fA,fB)[0]; autoA=S.cross(fA,fA)[0]; autoB=S.cross(fB,fB)[0]; full=S.cross(fF,fF)[0]; dd=S.cross(fD,fD)[0]
        T=S.theory(fA,fB,[Sl*pw**2])[0]; Tfsky=S.fsky_binned(Sl*pw**2)
        # Gaussian covariance of the half cross-spectrum from the measured half spectra
        cov=S.gaussian_covariance(fA,fB,fA,fB,S.spectrum_model(autoA,0.),S.spectrum_model(cross),S.spectrum_model(cross),S.spectrum_model(autoB,0.))
        band=bias_band(t.zmin,t.zmax); use=(S.ell_eff>=band[0])&(S.ell_eff<=band[1])&(T>0)
        b2,sb2,chi2,dof=fit_amplitude(cross,T,cov,use)
        b=float(np.sqrt(max(b2,1e-12))); sb=float(sb2/(2*b))
        # Shot noise of the full map = white level of the half difference. White pixel noise is LOCAL, so its level
        # is the pseudo-spectrum over <mask^2> (exact for white noise); the NaMaster decoupling of a spectrum that
        # is not band-limited attributes the coupling to ell > lmax to the band and overestimates the level by
        # 13-16 % for the fragmented DESI masks (2 % for BOSS). The decoupled value is kept as a diagnostic.
        pseudo=hp.anafast(diff,lmax=lmax,iter=0)/max(float(np.mean(mask.astype(float)**2)),1e-30)
        shot=float(np.mean(pseudo[100:lmax-40+1])); shot_decoupled=float(np.mean(dd[(S.ell_eff>=100)&(S.ell_eff<=lmax-40)]))
        fit={'b':b,'sigma_b':sb,'b2':float(b2),'sigma_b2':float(sb2),'chi2':chi2,'dof':dof,'band':[float(x) for x in band],'annuli_used':int(use.sum()),
             'estimator':'NaMaster cross-spectrum of two random halves (no shot noise), Gaussian covariance',
             'mask_transfer_band':(T/np.maximum(Tfsky,1e-30))[(S.ell_eff>=40)&(S.ell_eff<=300)].tolist(),
             'shot_from_half_difference':shot,'shot_from_half_difference_decoupled':shot_decoupled,'shot_model':meta['shot_s']}
        usable=np.isfinite(b2) and b2>3*sb2; b=fit['b'] if usable else t.bias
        fit['usable']=bool(usable); fit['bias_source']='auto-spectrum' if usable else 'fallback_table'
        maps.append(kmap/b); masks.append(mask); shots.append(shot/b**2); cats.append(cat); biases.append(b)
        uw.append(_unit_weights(cat,rnd,fp,cfg,nside,cref)/b)
        info['tracers'][t.label]={'n_objects':int(len(cat['ra'])),'sum_weights':float(cat['weight'].sum()),'fsky':fsky,'bias_fit':fit,'bias_used':float(b),
                                  'bias_table':t.bias,'shot_s_unit_bias':shot,'shot_s_uniform_model':meta['shot_s_uniform_model'],
                                  'L':S.ell_eff.tolist(),'ell_eff':S.ell_eff.tolist(),'cross_cl_binned':cross.tolist(),'cross_err':np.sqrt(np.diag(cov)).tolist(),
                                  'full_auto_binned':full.tolist(),'half_difference_binned':dd.tolist(),
                                  'theory_unit_bias_binned':T.tolist(),'theory_fsky_binned':Tfsky.tolist(),
                                  'data_random_ratio':meta['data_random_ratio'],'wall_s':time.perf_counter()-t0}
        hp.write_map(str(out/f'unitbias_{t.label}_nside{nside}.fits'),[kmap,mask.astype(float)],overwrite=True,dtype=np.float64)
        print(f"  {t.label}: {len(cat['ra'])} objects, fsky {fsky:.3f}, b = {fit['b']:.3f} +- {fit['sigma_b']:.3f} ({fit['bias_source']}, chi2 {chi2:.1f}/{dof}; "
              f"shot half-diff/model {shot/meta['shot_s']:.2f} (decoupled {shot_decoupled/meta['shot_s']:.2f}); transfer {np.round(fit['mask_transfer_band'][:3],2)}), {time.perf_counter()-t0:.0f} s",flush=True)
    # ---- noise matrix in convergence units: diagonal from the half difference, off-diagonal from shared objects
    k=len(maps); area=hp.nside2pixarea(nside); N=np.diag(shots); info['shared_objects']={}
    for i,j in combinations(range(k),2):
        joint=masks[i]&masks[j]
        if not joint.any(): continue
        ii,jj=shared_objects(cats[i],cats[j])
        nij=float(np.sum(uw[i][ii]*uw[j][jj])*area*area/max(float(joint.sum())*area,area))
        N[i,j]=N[j,i]=nij
        info['shared_objects'][f'{tracers[i].label}|{tracers[j].label}']={'n_shared':int(len(ii)),'fraction_of_first':float(len(ii)/max(len(cats[i]['ra']),1)),
                                                                            'fraction_of_second':float(len(ii)/max(len(cats[j]['ra']),1)),'noise_cross':nij,
                                                                            'noise_correlation':float(nij/np.sqrt(shots[i]*shots[j]))}
    # ---- Wiener combination per coverage class: C_kl = S pw^2 + N_kl on the tracers present, signal S pw
    union=np.zeros(len(masks[0]),bool)
    for m in masks: union|=m
    classes=coverage_classes(masks); comb=None; info['classes']={}; w_eff=np.zeros(lmax+1)
    for sub,cm in classes:
        idx=list(sub); W=np.zeros((lmax+1,len(idx)))
        for l in range(2,lmax+1):
            Cm=Sl[l]*pw[l]**2*np.ones((len(idx),len(idx)))+N[np.ix_(idx,idx)]
            W[l]=np.linalg.solve(Cm,np.full(len(idx),Sl[l]*pw[l]))
        frac=float(cm.sum()/union.sum()); w_eff+=frac*W.sum(axis=1)
        for a_,i in enumerate(idx):
            alm=hp.map2alm(maps[i]*cm,lmax=lmax,iter=0); part=hp.almxfl(alm,W[:,a_])
            comb=part if comb is None else comb+part
        info['classes']['+'.join(tracers[i].label for i in idx)]={'area_fraction':frac,'weights_at_L':{str(l):W[l].tolist() for l in (40,100,200,300,500)}}
    hp.write_alm(str(out/f"kappa_slice_{slice_['zmin']:g}_{slice_['zmax']:g}_alm.fits"),comb,overwrite=True)
    hp.write_map(str(out/f"mask_slice_{slice_['zmin']:g}_{slice_['zmax']:g}_nside{nside}.fits"),union.astype(float),overwrite=True,dtype=np.float64)
    info['shot_s']=shots; info['noise_matrix']=N.tolist(); info['fsky']=float(union.mean())
    info['effective_weight']=w_eff.tolist()      # area-weighted sum over tracers of the Wiener weights: <T kappa'> = w_eff C^{l c}
    info['weights_at_L']={str(l):[float(w_eff[l])] for l in (40,100,200,300)}
    return comb,union,info


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--nside',type=int,default=512); ap.add_argument('--lmax',type=int,default=1000)
    ap.add_argument('--tracers',nargs='*',default=['LRG','ELG','QSO','BGS','BOSS']); ap.add_argument('--out',type=Path,default=DATA/'lowz')
    ap.add_argument('--slices',nargs='*',type=float,default=None,help='zmin zmax pairs to restrict to')
    ap.add_argument('--zref',type=float,default=2.4,help='source-plane redshift of the forest (the weighted mean pixel redshift of the sample)')
    ap.add_argument('--tracer-zmax',type=float,default=1.75,help='drop tracer slices ending above this (iteration 11: 1.6, to keep a 300 Mpc/h buffer in front of a forest starting at z = 1.96)')
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True); cfg=Config(scale=1.,r_perp_min=3.,fit_rperp_min=3.).copy(chi_ref=float(chi_of_z(a.zref)))
    tracers=[t for t in TRACERS if t.name in a.tracers and t.zmax<=a.tracer_zmax+1e-6]; S=Spectra(a.lmax,width=int(ANNULUS))
    summary={'nside':a.nside,'lmax':a.lmax,'bin_width':int(ANNULUS),'spectra':'NaMaster','tracer_table':[t.label for t in tracers],'slices':[],
             'z_ref':a.zref,'chi_ref':float(cfg.chi_ref),'tracer_zmax':a.tracer_zmax}
    total=None; total_mask=None
    for s in slices_of(tuple(tracers)):
        if a.slices and not any(abs(s['zmin']-a.slices[i])<1e-6 and abs(s['zmax']-a.slices[i+1])<1e-6 for i in range(0,len(a.slices),2)): continue
        print(f"slice {s['zmin']}-{s['zmax']}: {[t.label for t in s['tracers']]}",flush=True)
        comb,mask,info=build_slice(s,s['tracers'],cfg,a.nside,a.lmax,a.out,S); summary['slices'].append(info)
        total=comb if total is None else total+comb; total_mask=mask if total_mask is None else (total_mask|mask)
        (a.out/'summary.json').write_text(json.dumps(summary,indent=1,default=float)+'\n')
    if total is not None:
        hp.write_alm(str(a.out/'kappa_combined_alm.fits'),total,overwrite=True)
        hp.write_map(str(a.out/f'mask_combined_nside{a.nside}.fits'),total_mask.astype(float),overwrite=True,dtype=np.float64)
    print('done',a.out)


if __name__=='__main__': main()
