"""Where does the lensing information come from? Fisher-information density of the pair estimator.

For a template with pair scalars d, the Gaussian independent-pair information is F = sum_pairs w_p w_q G^2 d^2
with G = chi_mid dxi/dr_perp. This script bins that sum (and the pair count and the mean G^2) in 1 Mpc/h cells of
(r_perp, r_par) for one saved scale-1 mock seed and its truth templates, and writes a JSON + figure.

Usage: python signal_profile.py --mock-root $LYALENSER_DATA/mocks/iteration4 --seed 0 [--output ../../report]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
from numba import njit, prange
HERE=Path(__file__).resolve().parent
for p in (HERE,HERE.parent):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
import run_mock_validation as v
from mock import sightlines_for_variant
from pairs import find_pairs, pair_geometry, _interp
from campaign4 import read_xi

@njit(parallel=True)
def _profile(pix_start, chi, weight, pa, pb, theta, dsq, rp_grid, rz_grid, xi, xirp, rpmax, rzmax, nbin):
    n=pa.size; nt=dsq.shape[0]
    counts=np.zeros((n,nbin,nbin)); gg=np.zeros((n,nbin,nbin)); ff=np.zeros((n,nt,nbin,nbin))
    rp0=rp_grid[0]; rz0=rz_grid[0]; drp=rp_grid[1]-rp_grid[0]; drz=rz_grid[1]-rz_grid[0]; nrp=rp_grid.size; nrz=rz_grid.size
    for ip in prange(n):
        a=pa[ip]; b=pb[ip]; qb=pix_start[b]; qend=pix_start[b+1]
        for p in range(pix_start[a],pix_start[a+1]):
            cp=np.float64(chi[p])
            while qb<qend and np.float64(chi[qb])<cp-rzmax: qb+=1
            q=qb
            while q<qend and np.float64(chi[q])<=cp+rzmax:
                cq=np.float64(chi[q]); rz=abs(cp-cq); cm=.5*(cp+cq); rp=cm*np.float64(theta[ip])
                if rp<rpmax and rz<rzmax:
                    xv,xg=_interp(rp,rz,rp0,drp,nrp,rz0,drz,nrz,xi,xirp)
                    G=cm*xg; ww=np.float64(weight[p])*np.float64(weight[q]); i=int(rp); j=int(rz)
                    counts[ip,i,j]+=1.; gg[ip,i,j]+=ww*G*G
                    for t in range(nt): ff[ip,t,i,j]+=ww*G*G*dsq[t,ip]
                q+=1
    return counts.sum(axis=0),gg.sum(axis=0),ff.sum(axis=0)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--mock-root',type=Path,required=True); ap.add_argument('--seed',type=int,default=0)
    ap.add_argument('--output',type=Path,default=v.ROOT/'report'); a=ap.parse_args()
    d=a.mock_root/f'sparse/{a.seed}'; m=v.load_mock(d/f'seed{a.seed:03d}.h5'); b=v.make_bundles(m)
    sl=sightlines_for_variant(m,0,False); xi=read_xi(d/f'fits{a.seed:03d}.h5','xi/A0_R0'); cfg=v.Config(scale=1)
    pa,pb,thx,thy,theta=find_pairs(sl,cfg.r_perp_max/sl.chi.min())
    science=[t for t in b['truth'] if getattr(t,'kind','') not in ('curl','junk')]
    if not science: raise RuntimeError('no science templates in the truth bundle: kinds '+str([getattr(t,'kind','') for t in b['truth']]))
    names=[f'{t.Lmin}-{t.Lmax}' for t in science]+['common science (sum of bands)']
    ds=[]
    for t in science:
        al=np.asarray(t.alpha,float); da=al[pa]-al[pb]; ds.append(thx*da[:,0]+thy*da[:,1])
    ds.append(np.sum(ds,axis=0)); dsq=np.asarray(ds)**2
    nb=int(cfg.r_perp_max)
    counts,gg,ff=_profile(sl.pix_start,sl.chi,sl.w,pa,pb,theta,dsq,xi.r_perp,xi.r_par,xi.xi.ravel(),xi.xi_rp.ravel(),cfg.r_perp_max,cfg.r_par_max,nb)
    total=ff[-1].sum(); frac_rp=ff[-1].sum(axis=1)/total; frac_rz=ff[-1].sum(axis=0)/total; cum=np.cumsum(frac_rp)
    out={'seed':a.seed,'sightlines':int(len(sl.ra)),'sightline_pairs':int(len(pa)),'pixel_pairs':float(counts.sum()),
         'bins_mpc':1.0,'templates':names,
         'information_fraction_vs_rperp':frac_rp.tolist(),'information_fraction_vs_rpar':frac_rz.tolist(),
         'cumulative_vs_rperp':cum.tolist(),'rperp_at_25_50_75_percent':[float(np.searchsorted(cum,x)+1) for x in (.25,.5,.75)],
         'per_band_fraction_vs_rperp':{n:(ff[i].sum(axis=1)/ff[i].sum()).tolist() for i,n in enumerate(names[:-1])},
         'pixel_pairs_vs_rperp':counts.sum(axis=1).tolist(),'mean_G2_vs_rperp':(gg.sum(axis=1)/np.maximum(counts.sum(axis=1),1)).tolist(),
         'information_2d':ff[-1].tolist(),'counts_2d':counts.tolist(),
         'note':'F = sum w w (chi_mid dxi/dr_perp)^2 d^2 per 1 Mpc/h cell; independent-pair Gaussian information, truth templates, response off, A_true=0 sample'}
    a.output.mkdir(parents=True,exist_ok=True); v.dump(a.output/'signal_profile.json',out)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    r=np.arange(nb)+.5; fig,ax=plt.subplots(1,3,figsize=(12,3.6))
    ax[0].bar(r,frac_rp,width=1,label='common science');
    for i,n in enumerate(names[:-1]): ax[0].plot(r,ff[i].sum(axis=1)/ff[i].sum(),label=f'L {n}')
    ax[0].set(xlabel='r_perp (Mpc/h)',ylabel='fraction of information per Mpc/h'); ax[0].legend(fontsize=8)
    ax[1].plot(r,counts.sum(axis=1)/counts.sum(),label='pixel pairs'); ax[1].plot(r,(gg.sum(axis=1)/np.maximum(counts.sum(axis=1),1))/np.max(gg.sum(axis=1)/np.maximum(counts.sum(axis=1),1)),label='mean G^2 (normalised)')
    ax[1].set(xlabel='r_perp (Mpc/h)',yscale='log'); ax[1].legend(fontsize=8)
    im=ax[2].imshow(ff[-1].T/total,origin='lower',extent=[0,nb,0,nb],aspect='auto'); ax[2].set(xlabel='r_perp',ylabel='r_par',title='information per cell'); fig.colorbar(im,ax=ax[2])
    fig.tight_layout(); (a.output/'figures').mkdir(exist_ok=True); fig.savefig(a.output/'figures/signal_profile.pdf'); plt.close(fig)
    print('information: 25/50/75% within r_perp =',out['rperp_at_25_50_75_percent'],'Mpc/h; r_par<10 share %.2f'%frac_rz[:10].sum())
    print('fraction per r_perp bin:',np.round(frac_rp,3))

if __name__=='__main__': main()
