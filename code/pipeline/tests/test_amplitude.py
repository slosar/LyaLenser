import numpy as np
from pairs import PairCatalogue
from templates import Template
from amplitude import amplitude,compress_score_per_sightline


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


def test_first_moment_formulas_and_known_ds():
    cat,t=synthetic(); g=.004; r=amplitude(cat,[t],g,np.array([0,1,2]))
    al=t.alpha.astype(float); d=cat.thx.astype(float)*(al[cat.a,0]-al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]-al[cat.b,1])
    s=cat.thx.astype(float)*(al[cat.a,0]+al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]+al[cat.b,1])
    x=cat.accum.sum(axis=2)
    q=np.sum((x[:,0]+g*x[:,1])*d+g*x[:,2]*s)
    mf=np.sum((x[:,8]+g*x[:,9])*d+g*x[:,10]*s)
    F=np.sum((x[:,3]+2*g*x[:,4]+g*g*x[:,5])*d*d+g*x[:,6]*d*s+g*g*x[:,7]*s*s)
    assert np.allclose([r.q[0],r.mf[0],r.F[0,0]],[q,mf,F],rtol=1e-14)
    assert np.isclose(r.A[0],(q-mf)/F)


def test_response_cross_template_symmetry():
    cat,t=synthetic(); u=Template(np.roll(t.alpha,1,axis=1),'other','curl')
    r=amplitude(cat,[t,u],.003,np.array([0,1,2]))
    assert np.allclose(r.F,r.F.T)


def test_per_sightline_compression_identity():
    cat,t=synthetic(); g=.003; U,score=compress_score_per_sightline(cat,t,g)
    x=cat.accum.sum(axis=2); h=x[:,0]-x[:,8]
    al=t.alpha.astype(float); d=cat.thx.astype(float)*(al[cat.a,0]-al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]-al[cat.b,1])
    s=cat.thx.astype(float)*(al[cat.a,0]+al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]+al[cat.b,1])
    direct=np.sum((h+g*(x[:,1]-x[:,9]))*d+g*(x[:,2]-x[:,10])*s)
    assert np.isclose(score,direct,rtol=1e-14)
