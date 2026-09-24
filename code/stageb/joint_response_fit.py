"""Joint response fit of all slice templates for the auto and the cross statistics, the normalised response-matrix
figure, and the auto x cross combination of the joint amplitudes (iteration 14).
Usage: python joint_response_fit.py --auto $LYALENSER_DATA/stageb/dr1_lowz_v7d --cross $LYALENSER_DATA/stageb/dr1_qso_v1d --lowz $LYALENSER_DATA/lowz_v4 --bands 40 200 400 600 800 1000 --tag v4
Writes report/stageb/joint_fit_<tag>.json and report/lowz/figures/response_matrix_<tag>.pdf.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
import h5py
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent; CODE=HERE.parent; ROOT=CODE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from campaign4 import campaign_config
from cosmo import chi as chi_of_z
from pairs import PairCatalogue, pair_midpoint_regions, ACCUMULATORS
from templates import sphere_band_templates, SCIENCE_BANDS
from joint_fit import build_joint, standard_fits, block_diagonal
from mock import load_sightlines
from qso_io import read_quasars
from xi_cross import Positions
from combine_auto_cross import combine


def load_cat(path,group='all'):
    with h5py.File(path) as f:
        g=f[group]; return PairCatalogue(g['a'][()],g['b'][()],g['thx'][()],g['thy'][()],g['theta'][()],g['accum'][()].astype(np.float64),g['npair'][()],{k:g.attrs[k] for k in g.attrs})


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--auto',type=Path,required=True); ap.add_argument('--cross',type=Path,default=None); ap.add_argument('--lowz',type=Path,required=True)
    ap.add_argument('--bands',type=float,nargs='+',default=None); ap.add_argument('--nside-alpha',type=int,default=2048); ap.add_argument('--nside-jk',type=int,default=8); ap.add_argument('--tag',default='joint')
    a=ap.parse_args(); t0=time.perf_counter()
    BANDS=tuple((int(a.bands[i]),int(a.bands[i+1])) for i in range(len(a.bands)-1)) if a.bands else SCIENCE_BANDS
    summary=json.loads((a.lowz/'summary.json').read_text()); slices=[f"slice_{s['zmin']:g}_{s['zmax']:g}" for s in summary['slices']]
    alms={n:hp.read_alm(str(a.lowz/f'kappa_{n}_alm.fits')) for n in slices}
    out={'bands':[list(b) for b in BANDS],'slices':slices,'statistics':{}}
    def run(kind,run,pos_ra,pos_dec,positions_for_regions,chi_ref):
        cfg=campaign_config(1.).copy(chi_ref=chi_ref); cat=load_cat(run/'catalogue.h5'); reg=pair_midpoint_regions(cat,positions_for_regions,a.nside_jk)
        print(f'[{kind}] {len(cat.a)} pairs, {len(np.unique(reg))} regions ({time.perf_counter()-t0:.0f} s)',flush=True)
        st={}
        for n in slices: st[n],_=sphere_band_templates(alms[n],pos_ra,pos_dec,nside=a.nside_alpha,science_bands=BANDS,source=n)
        print(f'[{kind}] templates evaluated ({time.perf_counter()-t0:.0f} s)',flush=True)
        jf=build_joint(cat,st,cfg.g1,reg); fits=standard_fits(jf); fits_bd=standard_fits(block_diagonal(jf))
        print(f"[{kind}] block-diagonal R: global A = {fits_bd['global']['A'][0]:.3f} +- {fits_bd['global']['jk_error'][0]:.3f}",flush=True)
        print(f"[{kind}] no curl: global A = {fits['no_curl_global']['A'][0]:.3f} +- {fits['no_curl_global']['jk_error'][0]:.3f}; no curl, no junk: {fits['no_curl_no_junk_global']['A'][0]:.3f} +- {fits['no_curl_no_junk_global']['jk_error'][0]:.3f}",flush=True)
        print(f"[{kind}] joint global A = {fits['global']['A'][0]:.3f} +- {fits['global']['jk_error'][0]:.3f} (F {fits['global']['sigma_F'][0]:.3f}); curl {fits['curl']['A'][0]:.3f} +- {fits['curl']['jk_error'][0]:.3f} ({time.perf_counter()-t0:.0f} s)",flush=True)
        for g,A,e in zip(fits['per_slice']['groups'],fits['per_slice']['A'],fits['per_slice']['jk_error']): print(f'[{kind}]   {g:16s} {A:7.3f} +- {e:.3f}',flush=True)
        for g,A,e in zip(fits['per_band']['groups'],fits['per_band']['A'],fits['per_band']['jk_error']): print(f'[{kind}]   {g:16s} {A:7.3f} +- {e:.3f}',flush=True)
        R=jf.normalised_response()
        rec={'n_components':len(jf.names),'names':jf.names,'kinds':jf.kinds,
             'global':{k:v for k,v in fits['global'].items() if k not in ('jk_samples','regions')},
             'per_slice':{k:v for k,v in fits['per_slice'].items() if k not in ('jk_samples','regions')},
             'per_band':{k:v for k,v in fits['per_band'].items() if k not in ('jk_samples','regions')},
             'curl':{k:v for k,v in fits['curl'].items() if k not in ('jk_samples','regions')},
             'block_diagonal_global':{k:v for k,v in fits_bd['global'].items() if k not in ('jk_samples','regions')},
             'no_curl_global':{k:v for k,v in fits['no_curl_global'].items() if k not in ('jk_samples','regions')},
             'no_curl_per_slice':{k:v for k,v in fits['no_curl_per_slice'].items() if k not in ('jk_samples','regions')},
             'no_curl_per_band':{k:v for k,v in fits['no_curl_per_band'].items() if k not in ('jk_samples','regions')},
             'no_curl_no_junk_global':{k:v for k,v in fits['no_curl_no_junk_global'].items() if k not in ('jk_samples','regions')},
             'response_normalised':R.tolist(),'max_offdiag_between_slices':float(np.max(np.abs(R[np.array([[n1.split(':')[0]!=n2.split(':')[0] for n2 in jf.names] for n1 in jf.names])])))}
        # separate-fit reference from the run's own products
        rr=json.loads((run/('dr1_lowz.json' if kind=='auto' else 'dr1_qso.json')).read_text()); rec['separate']={'combined':{k:rr['fits']['combined'][k] for k in ('A','jk_error','sigma_F')},'slices':{n:{k:rr['fits'][n][k] for k in ('A','jk_error','sigma_F')} for n in slices if n in rr['fits']},'bands':rr['fits']['combined']['bands'][:len(BANDS)],'band_errors':rr['fits']['combined']['band_errors'][:len(BANDS)]}
        rec['_fits_bd']=fits_bd
        return rec,fits,R
    sl=load_sightlines(a.auto/'sightlines.h5'); alog=json.loads((a.auto/'dr1_lowz.json').read_text()); chi_ref=float(alog['chi_ref'])
    rec,fa,Ra=run('auto',a.auto,sl.ra,sl.dec,sl,chi_ref); out['statistics']['auto']=rec
    if a.cross:
        clog=json.loads((a.cross/'dr1_qso.json').read_text()); qso=read_quasars(*clog['quasar_z'],verbose=False); pos=Positions(sl,qso)
        rec,fc,Rc=run('cross',a.cross,pos.ra,pos.dec,pos,chi_ref); out['statistics']['cross']=rec
        # combination of the joint amplitudes
        comb={'global':combine(fa['global']['A'][0],fa['global']['jk_samples'][:,0],fa['global']['regions'],fc['global']['A'][0],fc['global']['jk_samples'][:,0],fc['global']['regions'])}
        for key in ('no_curl_global','no_curl_no_junk_global'):
            comb[key]=combine(fa[key]['A'][0],fa[key]['jk_samples'][:,0],fa[key]['regions'],fc[key]['A'][0],fc[key]['jk_samples'][:,0],fc[key]['regions'])
        for key in ('per_slice','per_band','no_curl_per_slice','no_curl_per_band'):
            comb[key]={}
            for i,g in enumerate(fa[key]['groups']):
                comb[key][str(g)]=combine(fa[key]['A'][i],fa[key]['jk_samples'][:,i],fa[key]['regions'],fc[key]['A'][i],fc[key]['jk_samples'][:,i],fc[key]['regions'])
        fa_bd=out['statistics']['auto'].pop('_fits_bd'); fc_bd=out['statistics']['cross'].pop('_fits_bd')
        comb['block_diagonal_global']=combine(fa_bd['global']['A'][0],fa_bd['global']['jk_samples'][:,0],fa_bd['global']['regions'],fc_bd['global']['A'][0],fc_bd['global']['jk_samples'][:,0],fc_bd['global']['regions'])
        out['combination']=comb; c=comb['global']; print(f"[combined] joint: auto {c['auto'] if 'auto' in c else ''} -> A = {c['A']:.3f} +- {c['error']:.3f} (corr {c['correlation']:.2f}, weights {np.round(c['weights'],2)})",flush=True)
    out['statistics']['auto'].pop('_fits_bd',None); [out['statistics'][k].pop('_fits_bd',None) for k in out['statistics']]
    (ROOT/'report/stageb'/f'joint_fit_{a.tag}.json').write_text(json.dumps(out,indent=1,default=float)+'\n')
    plot_response(out,ROOT/'report/lowz/figures'/f'response_matrix_{a.tag}.pdf')
    print('saved',ROOT/'report/stageb'/f'joint_fit_{a.tag}.json')


def plot_response(out,path,vmax=0.3,only=None):
    """Normalised response matrices with thick lines between slices; the diagonal is one, the colour scale is
    +-vmax so that the off-diagonal structure is visible (it is a few per cent)."""
    mats=[(lab,np.asarray(out['statistics'][k]['response_normalised']),out['statistics'][k]['names']) for k,lab in (('auto','forest x forest'),('cross','quasar x forest')) if k in out['statistics'] and (only is None or k==only)]
    fig,axes=plt.subplots(1,len(mats),figsize=(6.2*len(mats),5.6)); axes=np.atleast_1d(axes)
    for ax,(title,R,names) in zip(axes,mats):
        n=len(names); Rp=R.copy(); np.fill_diagonal(Rp,np.nan); im=ax.imshow(Rp,vmin=-vmax,vmax=vmax,cmap='RdBu_r',origin='upper')
        sl_of=[nm.split(':')[0] for nm in names]; cuts=[i for i in range(1,n) if sl_of[i]!=sl_of[i-1]]
        for c in cuts: ax.axhline(c-.5,color='k',lw=2); ax.axvline(c-.5,color='k',lw=2)
        centres=[(np.mean([i for i,s in enumerate(sl_of) if s==sname])) for sname in dict.fromkeys(sl_of)]
        ax.set_xticks(centres); ax.set_xticklabels([s.replace('slice_','').replace('_','–') for s in dict.fromkeys(sl_of)],fontsize=8); ax.set_yticks(centres); ax.set_yticklabels([s.replace('slice_','').replace('_','–') for s in dict.fromkeys(sl_of)],fontsize=8)
        ax.set_title(f'{title}: $F_{{ij}}/\\sqrt{{F_{{ii}}F_{{jj}}}}$ ({n} components; diagonal masked)',fontsize=9); ax.grid(False)
    fig.colorbar(im,ax=axes.tolist(),shrink=.8); fig.savefig(path,bbox_inches='tight'); plt.close(fig)


if __name__=='__main__':
    if '--plot-only' in sys.argv:      # redraw from the saved product: python joint_response_fit.py --plot-only v4 [auto|cross]
        tag=sys.argv[sys.argv.index('--plot-only')+1]; only=sys.argv[-1] if sys.argv[-1] in ('auto','cross') else None
        out=json.loads((ROOT/'report/stageb'/f'joint_fit_{tag}.json').read_text()); suffix=f'_{only}' if only else ''
        plot_response(out,ROOT/'report/lowz/figures'/f'response_matrix{suffix}_{tag}.pdf',only=only); print('wrote',f'response_matrix{suffix}_{tag}.pdf')
    else: main()
