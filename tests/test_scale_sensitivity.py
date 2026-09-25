"""Check the diagnostic histogram against independently accumulated production pairs."""
import numpy as np
from lyalenser.config import Config, SightlineSet
from lyalenser.pairs import find_pairs, accumulate
from lyalenser.qso_io import QuasarSet
from lyalenser.xi_cross import find_cross_pairs, accumulate_cross
from lyalenser.xi_model import XiTable
from scripts.scale_sensitivity import _ff, _qf, response_total


def sample():
    sl=SightlineSet([1,2],[0.,.14],[0.,0.],[3.,3.],[0,2,4],
                   [3990.,4000.,4005.,4010.],np.ones(4),[1.,2.,3.,4.],[0]*4,{},[0,1,0,1])
    cfg=Config(chi_ref=3900.,g1=.001,r_perp_min=3.)
    grid=np.arange(0.,40.25,.25); shape=(len(grid),len(grid))
    tab=XiTable(grid,grid,np.ones(shape),np.ones(shape))
    return sl,cfg,tab


def test_auto_information_matches_production_response():
    sl,cfg,t=sample()
    cat=accumulate(sl,find_pairs(sl,.01),t,cfg)
    d=np.full(len(cat.a),.2); s=np.full(len(cat.a),.7)
    ff,count=_ff(sl.pix_start,sl.chi,sl.w,sl.slab,sl.region,cat.a,cat.b,cat.theta,
                 d,s,cfg.g1,cfg.chi_ref,t.r_perp,t.r_par,t.xi.ravel(),t.xi_rp.ravel(),
                 30.,30.,3.,*t.layers(),30,2)
    np.testing.assert_allclose(ff.sum(),response_total(cat,d,s,cfg.g1),rtol=1e-12)
    assert count.sum()==cat.npair.sum()==3


def test_cross_information_matches_production_response_and_boundary():
    sl,cfg,_=sample()
    q=QuasarSet(np.array([99]),np.array([.07]),np.array([0.]),np.array([2.4]),np.array([4030.],np.float32),{})
    rp=np.arange(0.,40.25,.25); rz=np.arange(-40.,40.25,.25)
    t=XiTable(rp,rz,np.ones((len(rp),len(rz))),np.ones((len(rp),len(rz))))
    cat=accumulate_cross(sl,q,find_cross_pairs(sl,q,.01),t,cfg)
    d=np.full(len(cat.a),.2); s=np.full(len(cat.a),.7)
    ff,count=_qf(sl.pix_start,sl.chi,sl.w,sl.slab,q.chi,cat.a,cat.b,cat.theta,sl.nq,
                 d,s,cfg.g1,cfg.chi_ref,t.r_perp,t.r_par,t.xi.ravel(),t.xi_rp.ravel(),
                 30.,30.,3.,*t.layers(),30,2)
    np.testing.assert_allclose(ff.sum(),response_total(cat,d,s,cfg.g1),rtol=1e-12)
    assert count.sum()==cat.npair.sum()==3
