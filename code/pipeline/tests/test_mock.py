import numpy as np
from scipy.ndimage import gaussian_filter
from mock import sample_lognormal_quasars


def test_lognormal_quasar_sampling_recovers_large_scale_bias():
    rng=np.random.default_rng(812)
    field=gaussian_filter(rng.normal(size=(64,64)),4,mode='wrap')
    field=.06*field/field.std()
    delta=np.repeat(field[:,:,None],4,axis=2).astype(np.float32)
    q=sample_lognormal_quasars(delta,2.,.5,3500.,100.,4000.,rng,b_q=3.5)
    counts=np.zeros((64,64),float)
    np.add.at(counts,(q['ix'],q['iy']),1)
    dq=counts/counts.mean()-1
    dq=gaussian_filter(dq,3,mode='wrap')
    dl=gaussian_filter(field,3,mode='wrap')
    bias=np.sum((dq-dq.mean())*(dl-dl.mean()))/np.sum((dl-dl.mean())**2)
    assert np.isclose(bias,3.5,rtol=.1),bias
