import numpy as np
from config import Config
from forest_power import ForestPower
from xi_model import xi_from_model,direct_variance


def test_kpar_normalisation_and_derivative():
    cfg=Config(xi_step=.5,xi_max=12)
    pf=ForestPower(model='kaiser',kp=5.)
    t=xi_from_model(pf,cfg,nk=360,kmax=18)
    independent=direct_variance(pf,nk=650,kmax=18)
    assert np.isclose(t.xi[0,0],independent,rtol=.01)
    fd=np.gradient(t.xi,cfg.xi_step,axis=0,edge_order=2)
    m=(t.r_perp>=1)&(t.r_perp<=8)
    scale=np.maximum(np.abs(t.xi_rp[m,:15]),1e-6*np.abs(t.xi[0,0]))
    assert np.nanmedian(np.abs(fd[m,:15]-t.xi_rp[m,:15])/scale)<.01


def test_grid_convergence_table_values():
    pf=ForestPower(model='kaiser',kp=4.)
    a=xi_from_model(pf,Config(xi_step=.5,xi_max=8),nk=320,kmax=15)
    b=xi_from_model(pf,Config(xi_step=.25,xi_max=8),nk=320,kmax=15)
    # Halving the real-space table step leaves coincident transform samples unchanged.
    assert np.allclose(a.xi,b.xi[::2,::2],rtol=1e-11,atol=1e-11)

