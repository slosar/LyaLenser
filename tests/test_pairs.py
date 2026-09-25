import numpy as np
import numba
from numba import set_num_threads
from lyalenser.config import Config,SightlineSet
from lyalenser.xi_model import XiTable
from lyalenser.pairs import find_pairs,accumulate,pair_geometry,ACCUMULATORS


def sample(nq=20,npix=50):
    rng=np.random.default_rng(4); c0=3900.
    ra=180+rng.uniform(-.15,.15,nq); dec=30+rng.uniform(-.15,.15,nq)
    chi=np.concatenate([c0+np.arange(npix)*.7+rng.uniform(0,.1) for _ in range(nq)]).astype('f4')
    return SightlineSet(np.arange(nq),ra,dec,np.full(nq,3,'f4'),np.arange(nq+1)*npix,chi,
                        rng.normal(size=nq*npix).astype('f4'),rng.uniform(.2,2,nq*npix).astype('f4'),
                        np.zeros(nq*npix,'i1'))


def table():
    g=np.arange(0,40.25,.25); xi=np.exp(-np.hypot(g[:,None],g[None,:])/12)
    return XiTable(g,g,xi,np.gradient(xi,.25,axis=0,edge_order=2))


def brute(sl,p,t,cfg):
    a,b,tx,ty,th=p; out=np.zeros((len(a),11,6)); count=np.zeros(len(a),int)
    for ip,(aa,bb) in enumerate(zip(a,b)):
      for i in range(sl.pix_start[aa],sl.pix_start[aa+1]):
       for j in range(sl.pix_start[bb],sl.pix_start[bb+1]):
        dc=float(sl.chi[i])-float(sl.chi[j]); rz=abs(dc); cm=.5*(float(sl.chi[i])+float(sl.chi[j])); rp=cm*th[ip]
        if rz>30 or rp>30: continue
        ir=0 if rp<10 else (1 if rp<20 else 2); iz=0 if rz<10 else 1; ib=2*ir+iz
        xv,xg=t.interp(rp,rz); G=cm*xg; dm=cm-cfg.chi_ref; ww=float(sl.w[i])*float(sl.w[j]); dd=float(sl.delta[i])*float(sl.delta[j])
        vals=(ww*dd*G,ww*dd*G*dm,ww*dd*G*dc/2,ww*G*G,ww*G*G*dm,ww*G*G*dm*dm,
              ww*G*G*dc,ww*G*G*dc*dc/4,ww*xv*G,ww*xv*G*dm,ww*xv*G*dc/2)
        out[ip,:,ib]+=vals; count[ip]+=1
    return out,count


def test_pair_kernel_vs_bruteforce():
    sl=sample(); cfg=Config(chi_ref=3900); p=find_pairs(sl,30/sl.chi.min()); t=table()
    cat=accumulate(sl,p,t,cfg); ref,n=brute(sl,p,t,cfg); keep=n>0
    assert np.allclose(cat.accum,ref[keep],rtol=1e-10,atol=1e-11)
    assert np.array_equal(cat.npair,n[keep])


def test_thread_determinism_bitwise():
    sl=sample(); cfg=Config(chi_ref=3900); p=find_pairs(sl,30/sl.chi.min()); t=table()
    set_num_threads(1); x=accumulate(sl,p,t,cfg).accum.copy()
    set_num_threads(numba.config.NUMBA_NUM_THREADS); y=accumulate(sl,p,t,cfg).accum.copy()
    assert np.array_equal(x,y)


def test_pair_direction_and_finite_difference_sign():
    ra=np.array([180.01,180.0]); dec=np.array([30.,30.]); a=np.array([0]); b=np.array([1])
    tx,ty,th=pair_geometry(ra,dec,a,b)
    assert tx[0]>0 and abs(ty[0])<1e-12
    chi=4000.; rp=chi*th[0]; kappa=1e-5
    xi=lambda r: np.exp(-r/10)
    finite=(xi((1-kappa)*rp)-xi(rp))/kappa
    alpha=np.array([[-kappa*th[0],0.],[0.,0.]])
    d=tx[0]*(alpha[0,0]-alpha[1,0]); response=chi*(-xi(rp)/10)*d/kappa
    assert finite>0 and np.isclose(response,finite,rtol=2e-5)
    assert -response<0  # reversing theta_ab flips the amplitude sign

