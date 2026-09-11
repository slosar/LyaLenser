import numpy as np
from pairs import PairCatalogue
from templates import Template
from amplitude import (_fit as amplitude,compress_score_per_sightline,
                       independent_response_prediction)


def synthetic():
    a=np.array([0,0,1],np.int32); b=np.array([1,2,2],np.int32)
    tx=np.array([1,0,2**-.5]); ty=np.array([0,1,2**-.5]); acc=np.zeros((3,11,6))
    vals=np.arange(1,19,dtype=float).reshape(3,6)
    acc[:,0]=vals; acc[:,1]=.2*vals; acc[:,2]=-.1*vals
    acc[:,3]=2+vals; acc[:,4]=.3*vals; acc[:,5]=.08*vals
    acc[:,6]=-.12*vals; acc[:,7]=.04*vals
    acc[:,8]=.5*vals; acc[:,9]=.07*vals; acc[:,10]=-.03*vals
    cat=PairCatalogue(a,b,tx,ty,np.ones(3),acc,np.ones(3,int))
    alpha=np.array([[.2,-.1],[-.1,.05],[.07,.2]])
    return cat,Template(alpha,'known','signal')


def with_junk(t):
    return [t,Template(np.roll(t.alpha,1,axis=1),"curl","curl"),
            Template(np.roll(t.alpha,1,axis=0),"junk","junk")]


def test_first_moment_formulas_and_known_ds():
    cat,t=synthetic(); g=.004; r=amplitude(cat,with_junk(t),g,np.array([0,1,2]))
    al=t.alpha.astype(float); d=cat.thx.astype(float)*(al[cat.a,0]-al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]-al[cat.b,1])
    s=cat.thx.astype(float)*(al[cat.a,0]+al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]+al[cat.b,1])
    x=cat.accum.sum(axis=2)
    q=np.sum((x[:,0]+g*x[:,1])*d+g*x[:,2]*s)
    mf=np.sum((x[:,8]+g*x[:,9])*d+g*x[:,10]*s)
    F=np.sum((x[:,3]+2*g*x[:,4]+g*g*x[:,5])*d*d+g*x[:,6]*d*s+g*g*x[:,7]*s*s)
    assert np.allclose([r.q[0],r.mf[0],r.F[0,0]],[q,mf,F],rtol=1e-14)
    assert np.allclose(r.F@r.A,r.q-r.mf)


def test_response_cross_template_symmetry():
    cat,t=synthetic(); u=Template(np.roll(t.alpha,1,axis=1),'other','curl')
    r=amplitude(cat,[t,u,Template(np.roll(t.alpha,1,axis=0),'junk','junk')],.003,np.array([0,1,2]))
    assert np.allclose(r.F,r.F.T)


def test_per_sightline_compression_identity():
    cat,t=synthetic(); g=.003; U,score=compress_score_per_sightline(cat,t,g)
    x=cat.accum.sum(axis=2); h=x[:,0]-x[:,8]
    al=t.alpha.astype(float); d=cat.thx.astype(float)*(al[cat.a,0]-al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]-al[cat.b,1])
    s=cat.thx.astype(float)*(al[cat.a,0]+al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]+al[cat.b,1])
    direct=np.sum((h+g*(x[:,1]-x[:,9]))*d+g*(x[:,2]-x[:,10])*s)
    assert np.isclose(score,direct,rtol=1e-14)


def test_response_matrix_with_signal_curl_and_junk():
    cat,t=synthetic()
    tc=Template(np.column_stack((-t.alpha[:,1],t.alpha[:,0])),'curl','curl')
    tj=Template(np.roll(t.alpha,1,axis=0),'junk','junk')
    r=amplitude(cat,[t,tc,tj],.002,np.array([0,1,2]))
    assert r.F.shape==(3,3) and np.allclose(r.F,r.F.T)
    assert np.linalg.matrix_rank(r.F)==3
    assert r.attrs['has_junk'] and r.attrs['has_curl']


def test_amplitude_refuses_missing_junk_band():
    import pytest
    cat,t=synthetic()
    with pytest.raises(ValueError,match='junk'):
        amplitude(cat,[t],0,np.array([0,1,2]))


def test_independent_response_prediction_synthetic():
    cat,density=synthetic()
    signal=Template(1.7*density.alpha,'signal','signal')
    tc=Template(np.column_stack((-density.alpha[:,1],density.alpha[:,0])),'curl','curl')
    tj=Template(np.roll(density.alpha,1,axis=0),'junk','junk')
    templates=[signal,tc,tj]; ratios=np.array([1.7,0.0,-.25])
    modulation=np.array([.2,-.1,.35])
    pred=independent_response_prediction(cat,templates,density,modulation,ratios,2.0,.003)
    expected_scores=pred['density_modulation_score']*ratios
    assert np.allclose(pred['scores'],expected_scores)
    assert np.allclose(pred['F']@pred['A'],expected_scores)
