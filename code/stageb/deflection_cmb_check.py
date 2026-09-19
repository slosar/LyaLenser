"""Validation of the deflection templates against CMB lensing (Stage B iteration 10, user request iv).

The pair estimator consumes the DEFLECTION alpha = grad phi of the Wiener-filtered tracer template (per slice and
combined), evaluated through the science-band windows. Here exactly that field is built on the sphere and
cross-correlated with the ACT DR6 and Planck PR4 convergence maps with NaMaster (spin-1 x spin-0): the E-mode of
the deflection carries the signal, the B-mode is a null. The prediction is not the kappa_lya x kappa_CMB spectrum
(the template estimates kappa_lya from tracers at z < 1.75 and is Wiener-suppressed where the tracers are noisy):
for a template T = sum_slices w_s(ell) * kappa_s_hat, <T kappa_CMB> = sum_s w_s(ell) C_ell^{(l_s, c)} with
C^{(l_s, c)} the Limber cross-spectrum of the slice's contribution to kappa_lya with kappa_CMB and w_s the
effective (area-weighted) Wiener weight of the slice (summary.json 'effective_weight'); the deflection E-mode
then reads C^{E kappa} = 2/sqrt(l(l+1)) F(l) <T kappa_CMB> for a band window F. The amplitude A_L of the measured
cross relative to this prediction is expected to be ONE if the template chain (biases, shot noise, Wiener
weights) is right; the ratio to the full kappa_lya x kappa_CMB spectrum is reported for reference.

Writes report/stageb/deflection_cmb_check.json.
Usage: python deflection_cmb_check.py [--lowz DIR] [--nside 512] [--lmax 1000] [--surveys ACT Planck]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from config import Config
from cosmo import chi as chi_of_z
from cross_spectrum import Z_CMB, kernel
from lowz import slice_spectra, ANNULUS
from templates import SCIENCE_BANDS, cosine_band
from nmt_spectra import Spectra, fit_amplitude, deflection_maps, gaussian_covariance_any
from cmb_maps import load_kappa, MASKED_ON_INPUT


def full_kappa_cross(cref,lmax):
    """kappa_lya x kappa_CMB (all redshifts) for reference."""
    from three_tracer import limber
    ccmb=float(chi_of_z(Z_CMB)); L=np.unique(np.r_[2.,np.linspace(2,lmax,400)])
    c=limber(L,lambda x: kernel(x,cref),lambda x: kernel(x,ccmb),1.,cref,nchi=600,to_recombination=False)
    return np.interp(np.arange(lmax+1),L,c)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--lowz',type=Path,default=DATA/'lowz_v2'); ap.add_argument('--nside',type=int,default=512)
    ap.add_argument('--lmax',type=int,default=1000); ap.add_argument('--out',type=Path,default=CODE.parent/'report/stageb/deflection_cmb_check.json')
    ap.add_argument('--surveys',nargs='*',default=['ACT','Planck']); ap.add_argument('--taper',type=float,default=10.)
    ap.add_argument('--bands',type=float,nargs='+',default=None,help='science band edges (default templates.SCIENCE_BANDS)')
    a=ap.parse_args(); cfg=Config(scale=1.,r_perp_min=3.,fit_rperp_min=3.); cref=cfg.chi_ref; lmax=a.lmax; ell=np.arange(lmax+1)
    summary=json.loads((a.lowz/'summary.json').read_text()); S=Spectra(lmax,width=int(ANNULUS)); pw=hp.pixwin(a.nside,lmax=lmax); cref=float(summary.get('chi_ref',cref))
    mask=hp.read_map(str(a.lowz/f'mask_combined_nside{a.nside}.fits'))
    # predictions per template (convergence units): sum_s w_eff,s C^{l_s c} pw, with the class-fraction average of the
    # Wiener weights taken INSIDE each CMB overlap (template_prediction.TemplatePrediction), not over the union
    from template_prediction import TemplatePrediction
    TP=TemplatePrediction(a.lowz,cref,lmax,a.nside)
    slices={}
    for s in summary['slices']:
        name=f"slice_{s['zmin']:g}_{s['zmax']:g}"; slices[name]=(a.lowz/f'kappa_{name}_alm.fits',hp.read_map(str(a.lowz/f'mask_{name}_nside{a.nside}.fits')))
    slices['combined']=(a.lowz/'kappa_combined_alm.fits',mask)
    full=full_kappa_cross(cref,lmax)
    BANDS=tuple((int(a.bands[i]),int(a.bands[i+1])) for i in range(len(a.bands)-1)) if a.bands else SCIENCE_BANDS
    windows={f'L{lo}_{hi}':cosine_band(ell,lo,hi,a.taper) for lo,hi in BANDS}
    windows['science']=sum(windows.values()); sl2=np.sqrt(np.maximum(ell*(ell+1.),1.)); grad=np.r_[0.,0.,2./sl2[2:]]
    out={'nside':a.nside,'lmax':lmax,'bands':{f'L{lo}_{hi}':[lo,hi] for lo,hi in BANDS},'surveys':{}}
    for survey in a.surveys:
        t0=time.perf_counter(); kmap,mk=load_kappa(survey,a.nside,lmax); fk=S.field(mk,[kmap],key=f'{survey}_mask',masked_on_input=MASKED_ON_INPUT[survey]); kk=S.cross(fk,fk)[0]
        res={'templates':{}}
        pred={name:TP.cross(tmask*mk,name) for name,(path,tmask) in slices.items()}
        res['effective_weight_on_overlap']={name:{str(l):float(TP.effective_weight(tmask*mk,name)[l]) for l in (40,100,200,300,500)} for name,(path,tmask) in slices.items() if name!='combined'}
        for name,(path,tmask) in slices.items():
            alm=hp.read_alm(str(path)); lm=hp.Alm.getlmax(len(alm)); r={'bands':{}}
            # spin-0 reference: the filtered kappa template itself against kappa_CMB
            fT=S.field(tmask,[hp.alm2map(hp.almxfl(alm,(np.arange(lm+1)<=lmax).astype(float)),a.nside,verbose=False)],key=f'{name}_mask')
            cx=S.cross(fT,fk)[0]; T=S.theory(fT,fk,[pred[name]])[0]; tt=S.cross(fT,fT)[0]
            cov=S.gaussian_covariance(fT,fk,fT,fk,S.spectrum_model(tt,0.),S.spectrum_model(cx),S.spectrum_model(cx),S.spectrum_model(kk,0.))
            use=(S.ell_eff>=BANDS[0][0])&(S.ell_eff<=BANDS[-1][1])&(T!=0); A,sA,chi2,dof=fit_amplitude(cx,T,cov,use)
            r['kappa_spin0']={'A':A,'sigma_A':sA,'chi2':chi2,'dof':dof,'L':S.ell_eff.tolist(),'cross':cx.tolist(),'cross_err':np.sqrt(np.diag(cov)).tolist(),
                              'prediction':T.tolist(),'prediction_full_kappa':S.theory(fT,fk,[full*pw])[0].tolist()}
            print(f"[{survey}] {name:16s} kappa x kappa_CMB ({BANDS[0][0]}-{BANDS[-1][1]}): A = {A:.3f} +- {sA:.3f}, chi2 {chi2:.1f}/{dof} ({time.perf_counter()-t0:.0f} s)",flush=True)
            # the deflection, per science band and for the science window
            for wname in list(windows):
                filt=windows[wname]; dth,dph=deflection_maps(alm,a.nside,lmax,filt=filt)
                fD=S.field(tmask,[dth,dph],spin=1,key=f'{name}_mask')
                cE,cB=S.cross(fD,fk); th=S.theory(fD,fk,[grad*filt*pred[name],np.zeros(lmax+1)]); TE=th[0]
                dd=S.cross(fD,fD)            # EE, EB, BE, BB
                cov4=gaussian_covariance_any(S,fD,fk,fD,fk,[S.spectrum_model(dd[0],0.),S.spectrum_model(dd[1]),S.spectrum_model(dd[2]),S.spectrum_model(dd[3],0.)],
                                             [S.spectrum_model(cE),S.spectrum_model(cB)],[S.spectrum_model(cE),S.spectrum_model(cB)],[S.spectrum_model(kk,0.)])
                covE=cov4[:,0,:,0]; covB=cov4[:,1,:,1]
                lo,hi=(BANDS[0][0],BANDS[-1][1]) if wname=='science' else tuple(int(x) for x in wname[1:].split('_'))
                use=(S.ell_eff>=lo-20)&(S.ell_eff<=hi+20)&(np.abs(TE)>0)&(np.abs(TE)>1e-3*np.abs(TE).max())
                A,sA,chi2,dof=fit_amplitude(cE,TE,covE,use)
                # B-mode: the deflection is a pure gradient, so its B-mode is mask leakage of the E-mode (deterministic,
                # not noise); its cross with kappa is quoted relative to the E-mode amplitude, with the E-mode
                # covariance as the scale (the Gaussian B covariance is ~0 and would give a meaningless error)
                Bn,sB,chi2B,dofB=fit_amplitude(cB,TE,covE,use)
                # the ratio to the full kappa_lya x kappa_CMB prediction (all redshifts, no Wiener suppression),
                # measured and expected with the same weighting
                TEf=S.theory(fD,fk,[grad*filt*full*pw,np.zeros(lmax+1)])[0]; Af,sAf,_,_=fit_amplitude(cE,TEf,covE,use); Ef,_,_,_=fit_amplitude(TE,TEf,covE,use)
                r['bands'][wname]={'A':A,'sigma_A':sA,'chi2':chi2,'dof':dof,'B_over_E':Bn/A if A!=0 else None,'B_amplitude':Bn,'B_amplitude_scale':sB,
                                   'A_vs_full_kappa_lya':Af,'sigma_A_vs_full':sAf,'expected_ratio_to_full':float(Ef),
                                   'L':S.ell_eff.tolist(),'E_cross':cE.tolist(),'E_err':np.sqrt(np.diag(covE)).tolist(),'B_cross':cB.tolist(),'B_err':np.sqrt(np.diag(covB)).tolist(),
                                   'prediction':TE.tolist(),'prediction_full_kappa_lya':TEf.tolist(),'used':use.tolist()}
                print(f"[{survey}] {name:16s} deflection {wname:9s}: A = {A:6.3f} +- {sA:.3f} (chi2 {chi2:.1f}/{dof}); B/E {Bn/A if A else float('nan'):6.3f}; vs full kappa_lya {Af:.3f} +- {sAf:.3f} (expected {Ef:.3f})",flush=True)
            res['templates'][name]=r
        res['wall_s']=time.perf_counter()-t0; out['surveys'][survey]=res
        a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(out,indent=1,default=float)+'\n')
    print('saved',a.out)


if __name__=='__main__': main()
