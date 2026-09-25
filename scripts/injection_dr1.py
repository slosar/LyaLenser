"""Injection tests for the auto and the cross statistics on the DR1 products (iteration 14): the combined template's
science deflection is injected by shifting every position (sightlines, and quasars for the cross) by -A alpha, the
pairs are rebuilt with the production selection and the amplitude refitted, for A = -0.5, -0.25, +0.25, +0.5.
Two modes per statistic: the noise-free expectation (delta products replaced by the fitted xi at the true
separation) and the actual injection into the real data; the recovered odd slope is reported for both.
Usage: python injection_dr1.py --auto $LYALENSER_DATA/stageb/dr1_lowz_v7d --cross $LYALENSER_DATA/stageb/dr1_qso_v1d --lowz $LYALENSER_DATA/lowz_v4 --bands 40 200 400 600 800 1000 --tag v4
Writes results/injection_<tag>.json.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
ROOT=Path(__file__).resolve().parents[1]
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
from lyalenser.config import production_config
from lyalenser.tables import read_xi
from lyalenser.pairs import find_pairs, accumulate, pair_midpoint_regions
from lyalenser.templates import sphere_band_templates, SCIENCE_BANDS
from lyalenser.amplitude import amplitude, curl_amplitude
from lyalenser.inject import shift_positions, paired_slopes, paired_jackknife
from lyalenser.xi_cross import find_cross_pairs, accumulate_cross, Positions
from lyalenser.desi_io import load_sightlines
from lyalenser.qso_io import read_quasars, QuasarSet
from lyalenser.amplitude import common_science
from lyalenser.tables import correction_for


def fit_A(cat,templates,cfg,reg):
    r=amplitude(cat,templates,cfg.g1,reg)
    s=common_science(r)
    n=sum(t.kind=='signal' for t in templates)
    return s['A'],curl_amplitude(r)[0],s['jk'],r.jk_samples[:,n:2*n].mean(axis=1),r.regions


def lensable_table(tab,basis_path,cfg):
    """Copy of the layered table with the same-wavelength term N(r_perp) 1[r_par < 1] subtracted from xi on every
    layer (it is not evolved with redshift; xi_rp already excludes it), rebuilt from the fit's coefficients with the
    production correction (`tables.correction_for`)."""
    from lyalenser.xi_fit import BASIS
    basis={'projected':{k:read_xi(basis_path,f'projected/{k}') for k in BASIS}}; corr=correction_for(basis,cfg)
    c=np.asarray(tab.meta['fit']['coefficients'],float); csw=np.r_[np.zeros(corr.n_s),c[corr.n_s:]]
    RP,RZ=np.meshgrid(tab.r_perp,tab.r_par,indexing='ij'); sw=corr.design(RP,RZ)@csw
    xi=np.asarray(tab.xi,float).copy(); xi-=sw[None,:,:] if xi.ndim==3 else sw
    from lyalenser.xi_model import XiTable
    return XiTable(tab.r_perp,tab.r_par,xi,tab.xi_rp,dict(tab.meta or {},note='same-wavelength term removed'),chi_nodes=getattr(tab,'chi_nodes',None))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--auto',type=Path,required=True); ap.add_argument('--cross',type=Path,default=None); ap.add_argument('--lowz',type=Path,required=True)
    ap.add_argument('--bands',type=float,nargs='+',default=None); ap.add_argument('--nside-alpha',type=int,default=2048); ap.add_argument('--nside-jk',type=int,default=8)
    ap.add_argument('--amplitudes',type=float,nargs='+',default=[-.5,-.25,.25,.5]); ap.add_argument('--tag',default='inj'); ap.add_argument('--modes',nargs='+',default=['expectation','real'])
    ap.add_argument('--statistics',nargs='+',choices=['auto','cross'],default=['auto','cross'],help='statistics to run; cross also requires --cross')
    ap.add_argument('--remove-same-wavelength',action='store_true',help='expectation with the same-wavelength term N(r_perp) removed from the injected correlation (the kernel omits it)'); ap.add_argument('--basis',type=Path,default=DATA/'stageb/basis_dr1_ab_z196.h5')
    a=ap.parse_args(); t0=time.perf_counter()
    BANDS=tuple((int(a.bands[i]),int(a.bands[i+1])) for i in range(len(a.bands)-1)) if a.bands else SCIENCE_BANDS
    alm=hp.read_alm(str(a.lowz/'kappa_combined_alm.fits')); out={'bands':[list(b) for b in BANDS],'amplitudes':a.amplitudes,'statistics':{}}
    out['fit']='combined template, five science amplitudes constrained equal; five curl and one junk nuisance'
    out['region_labels_preserved']=True
    def checkpoint():
        (ROOT/'results'/f'injection_{a.tag}.json').write_text(json.dumps(out,indent=1,default=float)+'\n')
    alog=json.loads((a.auto/'dr1_lowz.json').read_text()); chi_ref=float(alog['chi_ref']); sl=load_sightlines(a.auto/'sightlines.h5')
    # ---- auto
    cfg=production_config().copy(chi_ref=chi_ref,slabs=((alog['zmin'],alog['zmax']),),r_perp_max=float(alog['config'].get('r_perp_max',30.))); tab=read_xi(a.auto/'xi.h5','xi')
    if a.remove_same_wavelength: tab=lensable_table(tab,a.basis,cfg); out['note']='same-wavelength term removed from the injected correlation (expectation mode)'
    templates,_=sphere_band_templates(alm,sl.ra,sl.dec,nside=a.nside_alpha,science_bands=BANDS,source='combined'); alpha_inj=sum(t.alpha for t in templates if t.kind=='signal')
    rec={}
    for mode in a.modes if 'auto' in a.statistics else []:
        vals=[]; curls=[]; samples=[]; curl_samples=[]; regions=[]; counts=[]
        for A in a.amplitudes:
            shifted=shift_positions(sl,alpha_inj,A); ps=find_pairs(shifted,cfg.r_perp_max/max(float(shifted.chi.min()),1))
            cat=accumulate(shifted,ps,tab,cfg,true_positions=np.column_stack((sl.ra,sl.dec)) if mode=='expectation' else None); reg=pair_midpoint_regions(cat,shifted,a.nside_jk)
            Ah,cu,jk,cjk,regs=fit_A(cat,templates,cfg,reg)
            vals.append(Ah); curls.append(cu); samples.append(jk); curl_samples.append(cjk); regions.append(regs); counts.append(int(cat.npair.sum()))
            print(f'[auto {mode}] A_inj {A:+.2f}: A_hat {Ah:.4f}, curl {cu:.3f} ({time.perf_counter()-t0:.0f} s)',flush=True)
        slope,per=paired_slopes(np.asarray(a.amplitudes),np.asarray(vals)); rec[mode]={'A_hat':vals,'curl':curls,'paired_slope':slope,'paired_slopes_by_amplitude':per}
        rec[mode].update(paired_jackknife(a.amplitudes,vals,samples,regions))
        rec[mode]['curl_slope']=paired_slopes(a.amplitudes,curls)[0]
        rec[mode]['curl_slope_jk_error']=paired_jackknife(a.amplitudes,curls,curl_samples,regions)['paired_slope_jk_error']
        rec[mode]['accepted_pixel_pairs']=counts
        out['statistics']['auto']=rec
        checkpoint()
        print(f'[auto {mode}] odd slope {slope:.4f} {per}',flush=True)
    # ---- cross
    if a.cross and 'cross' in a.statistics:
        clog=json.loads((a.cross/'dr1_qso.json').read_text()); qso=read_quasars(*clog['quasar_z'],verbose=False); pos=Positions(sl,qso); tabc=read_xi(a.cross/'xi_qf.h5','xi')
        cfgc=production_config().copy(chi_ref=chi_ref,slabs=((alog['zmin'],alog['zmax']),),r_perp_max=float(clog['config'].get('r_perp_max',30.)))
        tpl,_=sphere_band_templates(alm,pos.ra,pos.dec,nside=a.nside_alpha,science_bands=BANDS,source='combined'); alpha_all=sum(t.alpha for t in tpl if t.kind=='signal')
        rec={}
        for mode in a.modes:
            vals=[]; curls=[]; samples=[]; curl_samples=[]; regions=[]; counts=[]
            for A in a.amplitudes:
                sh=shift_positions(sl,alpha_all[:sl.nq],A); dec_rad=np.deg2rad(qso.dec)
                qsh=QuasarSet(qso.qid,qso.ra-np.rad2deg(A*alpha_all[sl.nq:,0]/np.maximum(np.cos(dec_rad),1e-8)),qso.dec-np.rad2deg(A*alpha_all[sl.nq:,1]),qso.z,qso.chi,qso.attrs)
                ps=find_cross_pairs(sh,qsh,cfgc.r_perp_max/max(1.,float(min(sh.chi.min(),qsh.chi.min()))))
                cat=accumulate_cross(sh,qsh,ps,tabc,cfgc,true_positions=np.column_stack((pos.ra,pos.dec)) if mode=='expectation' else None); reg=pair_midpoint_regions(cat,Positions(sh,qsh),a.nside_jk)
                Ah,cu,jk,cjk,regs=fit_A(cat,tpl,cfgc,reg)
                vals.append(Ah); curls.append(cu); samples.append(jk); curl_samples.append(cjk); regions.append(regs); counts.append(int(cat.npair.sum()))
                print(f'[cross {mode}] A_inj {A:+.2f}: A_hat {Ah:.4f}, curl {cu:.3f} ({time.perf_counter()-t0:.0f} s)',flush=True)
            slope,per=paired_slopes(np.asarray(a.amplitudes),np.asarray(vals)); rec[mode]={'A_hat':vals,'curl':curls,'paired_slope':slope,'paired_slopes_by_amplitude':per}
            rec[mode].update(paired_jackknife(a.amplitudes,vals,samples,regions))
            rec[mode]['curl_slope']=paired_slopes(a.amplitudes,curls)[0]
            rec[mode]['curl_slope_jk_error']=paired_jackknife(a.amplitudes,curls,curl_samples,regions)['paired_slope_jk_error']
            rec[mode]['accepted_pixel_pairs']=counts
            out['statistics']['cross']=rec
            checkpoint()
            print(f'[cross {mode}] odd slope {slope:.4f} {per}',flush=True)
        out['statistics']['cross']=rec
    checkpoint(); print('saved')


if __name__=='__main__': main()
