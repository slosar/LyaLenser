import numpy as np
from lyalenser.config import SightlineSet
from lyalenser.inject import shift_positions


def test_injection_shift_is_minus_alpha():
    sl=SightlineSet([1],[10.],[30.],[3.],[0,1],[4000.],[0.],[1.],[0])
    alpha=np.array([[1e-4,2e-4]])
    out=shift_positions(sl,alpha,2)
    assert out.dec[0]<sl.dec[0] and out.ra[0]<sl.ra[0]
    assert np.isclose(np.deg2rad(out.dec[0]-sl.dec[0]),-4e-4)

