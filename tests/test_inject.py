import numpy as np
from lyalenser.config import SightlineSet
from lyalenser.inject import shift_positions


def test_injection_shift_is_minus_alpha():
    sl=SightlineSet([1],[10.],[30.],[3.],[0,1],[4000.],[0.],[1.],[0])
    alpha=np.array([[1e-4,2e-4]])
    out=shift_positions(sl,alpha,2)
    assert out.dec[0]<sl.dec[0] and out.ra[0]<sl.ra[0]
    assert np.isclose(np.deg2rad(out.dec[0]-sl.dec[0]),-4e-4)


def test_zero_shift_preserves_bb_exclusion():
    from lyalenser.config import Config
    from lyalenser.pairs import find_pairs, accumulate
    from lyalenser.xi_model import XiTable
    sl=SightlineSet([1,2],[0.,.14],[0.,0.],[3.,3.],[0,1,2],
                   [4000.,4000.],[1.,1.],[1.,1.],[0,0],{},[1,1])
    shifted=shift_positions(sl,np.zeros((2,2)),0.)
    np.testing.assert_array_equal(shifted.region,sl.region)
    assert not np.shares_memory(shifted.region,sl.region)
    cfg=Config(chi_ref=4000.,g1=0.,r_perp_min=3.)
    grid=np.arange(41.); tab=XiTable(grid,grid,np.ones((41,41)),np.ones((41,41)))
    for sample in (sl,shifted):
        assert accumulate(sample,find_pairs(sample,.01),tab,cfg).npair.sum()==0


def test_paired_jackknife_cancels_common_noise():
    from lyalenser.inject import paired_jackknife
    regions=[np.arange(4),np.arange(4)]
    noise=np.array([.2,-.1,.4,-.5])
    result=paired_jackknife([-1.,1.],[-1.,1.],[-1+noise,1+noise],regions)
    assert result['paired_slope_jk_error'] < 1e-14
