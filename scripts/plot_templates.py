"""Figure: the Wiener-filtered convergence templates per slice and combined (iteration 14, paper figure).
Reads $LYALENSER_DATA/<lowz>/kappa_*_alm.fits and the masks; shows each map smoothed to the science window
(40 <= L <= 1000 cosine-tapered, then 20 arcmin Gaussian for display) in a Mollweide projection, equatorial.
Writes Paper/figures/templates_maps.pdf and report/figures/templates_maps.pdf.
Usage: python plot_templates.py --lowz $LYALENSER_DATA/lowz_v4
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import healpy as hp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA
from lyalenser.templates import cosine_band


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--lowz',type=Path,default=DATA/'lowz_v4'); ap.add_argument('--nside',type=int,default=256); ap.add_argument('--fwhm',type=float,default=30.); ap.add_argument('--lmax-show',type=int,default=400); ap.add_argument('--combined-only',action='store_true',help='one panel, the combined template only (paper Figure 4) -> templates_map_combined.pdf')
    a=ap.parse_args(); summary=json.loads((a.lowz/'summary.json').read_text()); slices=[f"slice_{s['zmin']:g}_{s['zmax']:g}" for s in summary['slices']]
    names=(['combined'] if a.combined_only else slices+['combined']); n=len(names); fig=plt.figure(figsize=(7,4.2) if a.combined_only else (12,2.9*((n+1)//2)))
    lmax_alm=hp.Alm.getlmax(len(hp.read_alm(str(a.lowz/'kappa_combined_alm.fits')))); ell=np.arange(lmax_alm+1); win=cosine_band(ell,40,a.lmax_show,10)
    for i,name in enumerate(names):
        alm=hp.read_alm(str(a.lowz/f'kappa_{name}_alm.fits')); alm=hp.almxfl(alm,win)
        m=hp.alm2map(alm,a.nside,fwhm=np.deg2rad(a.fwhm/60),verbose=False); mask=hp.ud_grade(hp.read_map(str(a.lowz/f"mask_{name}_nside{summary['nside']}.fits")),a.nside)>0.5
        m=np.where(mask,m,hp.UNSEEN); lim=np.percentile(np.abs(m[mask]),99)
        ttl=(name.replace('slice_','tracers at ').replace('_',' < z < ') if name!='combined' else 'combined template (sum of the slices)')
        hp.mollview(m,sub=(None if a.combined_only else ((n+1)//2,2,i+1)),title=('' if a.combined_only else ttl),min=-lim,max=lim,cmap='RdBu_r',cbar=False,notext=True,margins=(0.01,0.02,0.01,0.02),fig=fig.number)
    fig.text(0.5,0.02,r'$\hat\kappa_{\rm Ly\alpha}$, Wiener-filtered, $40\leq \ell\leq%d$, %g arcmin smoothing; colour scale $\pm$ the 99th percentile%s'%(a.lmax_show,a.fwhm,'' if a.combined_only else ' of each map'),ha='center',fontsize=9)
    base='templates_map_combined.pdf' if a.combined_only else 'templates_maps.pdf'
    for f in (ROOT/'Paper/figures'/base,ROOT/'report/figures'/base):
        f.parent.mkdir(parents=True,exist_ok=True); fig.savefig(f,bbox_inches='tight',dpi=150)
    print('wrote figures'); plt.close(fig)


if __name__=='__main__': main()
