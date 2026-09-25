"""Iteration-6 diagnosis (2026-09-14): the forest's own long-wavelength component in the pair products.

The scale-1 campaign (N = 400) shows a template-correlated null with the response OFF (varying-template mean
field: matched 10.3 +- 1.6, cmb 0.33 +- 0.14, deprojected -0.41 +- 0.18, all in A; the fixed-template rows are
zero) and a response 11-18 % above the P (D C D - C) P^T prediction (paired A0_R1 - A0_R0: matched 28.7 +- 0.7
vs 23.7 predicted, cmb 1.88 +- 0.06 vs 1.63). Hypothesis (NOTES.md, iteration 6): the pair product
delta_p delta_q contains m_p m_q, the product of the forest's long-mode components, which is coherent with the
template field; the mean field subtracts only xi_pq (a realisation average) and the prediction conditions on
delta_L while treating delta_F as independent of it. For a lognormal template <m_p m_q T> = b_q^2 xi_ps xi_qs
exactly, and with the modulation on the Gaussian four-point terms <delta_L m_p><m_q T> appear.

Test: m_p = c P[delta_L(pixel)] with P the per-forest continuum projection and c the weighted regression of the
A0_R0 forest on P[delta_L] (one scalar per seed, fixed for every variant, so the lensing and the response are not
refitted); delta' = delta - m. The kernel (xi_rp of the frozen model-fit table of the standard sample) is kept;
only the mean field changes, by the MEASURED difference of the 1 Mpc/h correlation of delta and delta'
(same pairs, same weights: the pixel-noise part cancels in the difference), turned into a fine table by the
campaign's even-reflection spline. (A refit of the two-parameter Kaiser table to delta' is not possible: the
long modes carry a large part of xi at 10-30 Mpc/h, and the refit runs away; first attempt, smoke seed.)
The m-only field (m_p m_q against its own measured xi, same kernel) isolates <m m T>. Every variant of the saved
seed is refitted with the seed's own templates. Reads campaign products, writes its own JSON per seed. Bypasses
the provenance check on purpose: this file changes the code fingerprint, and the campaign is complete.

Usage:
  python longmode_diagnosis.py --mock-root $LYALENSER_DATA/mocks/iteration6 --scale 1 --seed 1000 --output DIR
  python longmode_diagnosis.py --summarize DIR
"""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from config import SightlineSet
from mock import project_continuum,sightlines_for_variant
from campaign4 import campaign_config,read_xi,BASIS
from xi_model import XiTable,xi_from_data,xi_from_counts
import run_mock_validation as v

VARIANTS=((0,False),(0,True),(1,True),(1,False))
NAMES=('truth','cmb','matched','deprojected')

def load_basis(root):
    """campaign_basis without the provenance check."""
    import h5py
    path=root/'basis/basis.h5'
    out={'raw':{k:read_xi(path,f'raw/{k}') for k in BASIS},'projected':{k:read_xi(path,f'projected/{k}') for k in BASIS}}
    with h5py.File(path) as f: out['coarse']={k:f[f'coarse/{k}'][()] for k in BASIS}; out['cpix']=f['cpix'][()]
    return out

def projected_long(m):
    s=m.sightlines; nq=len(s.pix_start)-1; npc=len(s.chi)//nq
    if nq*npc!=len(s.chi): raise ValueError('forests of unequal length')
    dL=np.asarray(m.truth['delta_L_pixel'],float).reshape(nq,npc)
    return project_continuum(dL,np.asarray(s.chi,float).reshape(nq,npc),np.asarray(s.w,float).reshape(nq,npc)).ravel()

def with_delta(sl,delta,**attrs):
    return SightlineSet(sl.qid,sl.ra,sl.dec,sl.zq,sl.pix_start,sl.chi,np.asarray(delta,np.float32),sl.w,sl.slab,{**sl.attrs,**attrs})

def coarse_summary(num,den,ref_num=None):
    """Coarse xi along r_par < 1 (first r_par bin) at r_perp bins 3-30, and its ratio to a reference."""
    raw=np.divide(num,den,out=np.zeros_like(num),where=den>0); row=raw[3:30,0]
    out={'r_perp_bins':'3..29 (+0.5)','r_par_bin':'0-1','xi':row.tolist()}
    if ref_num is not None:
        ref=np.divide(ref_num,den,out=np.zeros_like(num),where=den>0)[3:30,0]
        out['ratio_to_reference']=np.divide(row,ref,out=np.zeros_like(row),where=ref!=0).tolist()
    return out

def run_seed(seed,root,scale,output):
    t=time.perf_counter(); cfg=campaign_config(scale); basis=load_basis(root)
    directory=root/f'sparse/{seed}'
    m=v.load_mock(directory/f'seed{seed:03d}.h5')
    standard=json.loads((directory/f'diagnostics{seed:03d}.json').read_text())
    b=v.make_bundles(m)
    pairs=v.find_pairs(m.sightlines,cfg.r_perp_max/m.sightlines.chi.min())
    P=projected_long(m); w=np.asarray(m.sightlines.w,float)
    d0=np.asarray(m.truth['delta_A0_response0'],float)
    c=float(np.sum(w*d0*P)/np.sum(w*P*P)); mfield=c*P
    fraction=float(np.sum(w*mfield**2)/np.sum(w*d0*d0))
    out={'seed':seed,'c':c,'variance_fraction_removed':fraction,'rms_long':float(np.sqrt(np.average(mfield**2,weights=w))),
         'rms_forest_A0_R0':float(np.sqrt(np.average(d0**2,weights=w))),
         'standard':{},'subtracted':{},'prediction_standard':{},'xi_fit_standard':{},'xi_change':{},'m_only':{}}
    # m-only field: its own measured xi for the mean field, the standard A0_R0 kernel.
    sl0=sightlines_for_variant(m,0,False); ft0=v.table_for(sl0,cfg,basis)
    slm=with_delta(sl0,mfield,long_mode_only=c); tm=xi_from_data(slm,cfg)
    table_m=XiTable(ft0.table.r_perp,ft0.table.r_par,tm.xi,ft0.table.xi_rp,{'provider':'m-only measured xi, standard kernel'})
    catm=v.cat_for(slm,table_m,cfg,pairs)
    out['m_only']={n:v.fit_bundle(catm,b[n],cfg.g1,slm)['A'] for n in NAMES}
    out['xi_m_only']=coarse_summary(*tm.counts,ref_num=ft0.table.counts[0])
    print('m-only',out['m_only'],flush=True)
    for A,resp in VARIANTS:
        key=f'A{A:g}_R{int(resp)}'
        if key not in standard['fits']: continue
        out['standard'][key]={n:standard['fits'][key][n]['A'] for n in NAMES}
        out['xi_fit_standard'][key]=standard['xi_fit'][key]
        sl=sightlines_for_variant(m,A,resp); ft=ft0 if (A,resp)==(0,False) else v.table_for(sl,cfg,basis)
        num,den=ft.table.counts
        sl2=with_delta(sl,np.asarray(sl.delta,float)-mfield,long_mode_subtracted=c)
        num2,den2=xi_from_data(sl2,cfg).counts
        if not np.allclose(den,den2,rtol=1e-6,atol=0): raise RuntimeError('pair weights differ between delta and delta-m')
        dxi=xi_from_counts(num-num2,den,cfg)
        out['xi_change'][key]=coarse_summary(num-num2,den,ref_num=num)
        table2=XiTable(ft.table.r_perp,ft.table.r_par,ft.table.xi-dxi.xi,ft.table.xi_rp,
                       {'provider':'standard model-fit xi minus measured (xi_delta - xi_delta-m), standard kernel'})
        cat=v.cat_for(sl2,table2,cfg,pairs)
        out['subtracted'][key]={n:v.fit_bundle(cat,b[n],cfg.g1,sl2)['A'] for n in NAMES}
        if A==0 and resp:
            out['prediction_standard']={n:standard['predictions'][n]['A'] for n in ('cmb','matched','deprojected')}
        print(key,'standard',out['standard'][key],'subtracted',out['subtracted'][key],flush=True)
    out['wall_s']=time.perf_counter()-t
    output.mkdir(parents=True,exist_ok=True)
    (output/f'longmode{seed}.json').write_text(json.dumps(out,indent=1,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x))+'\n')
    print('DONE',seed,'c',c,'fraction',fraction,'wall',out['wall_s'],flush=True)
    return out

def stats(x):
    x=np.asarray(x,float); return {'n':int(len(x)),'mean':float(x.mean()),'sem':float(x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else float('nan')}

def summarize(output):
    files=sorted(Path(output).glob('longmode*.json')); R=[json.loads(f.read_text()) for f in files]
    print(f'{len(R)} seeds; c mean {np.mean([r["c"] for r in R]):.4f} +- {np.std([r["c"] for r in R]):.4f}; '
          f'variance fraction removed {np.mean([r["variance_fraction_removed"] for r in R]):.4f}; '
          f'rms m {np.mean([r["rms_long"] for r in R]):.4f}, rms forest {np.mean([r["rms_forest_A0_R0"] for r in R]):.4f}')
    S={}
    def line(label,std,sub):
        s1,s2,sd=stats(std),stats(sub),stats(np.asarray(std)-np.asarray(sub))
        S[label]={'standard':s1,'subtracted':s2,'paired_difference':sd}
        print(f'{label:52s} standard {s1["mean"]:8.3f} +- {s1["sem"]:6.3f} | subtracted {s2["mean"]:8.3f} +- {s2["sem"]:6.3f} | std-sub {sd["mean"]:8.3f} +- {sd["sem"]:6.3f}')
    for n in NAMES:
        g=lambda k,which:np.array([r[which][k][n] for r in R if k in r[which]])
        line(f'{n} A0_R0 (response off, own templates)',g('A0_R0','standard'),g('A0_R0','subtracted'))
        mo=stats([r['m_only'][n] for r in R]); S[f'{n} m-only']=mo
        print(f'{n+" m-only field":52s} {mo["mean"]:8.3f} +- {mo["sem"]:6.3f}')
        line(f'{n} A0_R1 - A0_R0 (response)',g('A0_R1','standard')-g('A0_R0','standard'),g('A0_R1','subtracted')-g('A0_R0','subtracted'))
        if n!='truth':
            p1=np.array([r['prediction_standard'][n] for r in R])
            line(f'{n} response minus prediction',g('A0_R1','standard')-g('A0_R0','standard')-p1,g('A0_R1','subtracted')-g('A0_R0','subtracted')-p1)
        line(f'{n} A1_R1 - A0_R1 (lensing, response on)',g('A1_R1','standard')-g('A0_R1','standard'),g('A1_R1','subtracted')-g('A0_R1','subtracted'))
        if all('A1_R0' in r['standard'] for r in R):
            line(f'{n} A1_R0 - A0_R0 (lensing, response off)',g('A1_R0','standard')-g('A0_R0','standard'),g('A1_R0','subtracted')-g('A0_R0','subtracted'))
        line(f'{n} A1_R1 (combined recovery)',g('A1_R1','standard'),g('A1_R1','subtracted'))
    S['seeds']=[r['seed'] for r in R]; S['c']=stats([r['c'] for r in R]); S['variance_fraction_removed']=stats([r['variance_fraction_removed'] for r in R])
    ratio=np.array([r['xi_change']['A0_R0']['ratio_to_reference'] for r in R]); S['xi_change_ratio_A0_R0_rpar0']={'r_perp_bins':'3..29','mean':ratio.mean(0).tolist()}
    print('(xi_delta - xi_delta-m)/xi_delta at r_par<1, r_perp 3..29:',np.round(ratio.mean(0),3))
    ratio=np.array([r['xi_m_only']['ratio_to_reference'] for r in R]); S['xi_m_only_ratio_rpar0']={'r_perp_bins':'3..29','mean':ratio.mean(0).tolist()}
    print('xi_mm/xi_delta at r_par<1, r_perp 3..29:',np.round(ratio.mean(0),3))
    (Path(output)/'summary.json').write_text(json.dumps(S,indent=1)+'\n')
    return S

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--mock-root',type=Path); ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--seed',type=int,nargs='*',default=[]); ap.add_argument('--output',type=Path)
    ap.add_argument('--summarize',type=Path)
    a=ap.parse_args()
    if a.summarize: summarize(a.summarize)
    else:
        for s in a.seed: run_seed(s,a.mock_root,a.scale,a.output)
