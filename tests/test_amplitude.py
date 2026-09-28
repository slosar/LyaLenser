import numpy as np
from lyalenser.pairs import PairCatalogue
from lyalenser.templates import Template
from lyalenser.amplitude import (_fit as amplitude,compress_score_per_sightline,
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
    cat,t=synthetic(); g=.003; U,Ud,score=compress_score_per_sightline(cat,t,g)
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


def pair_scalar_pairs(cat,al):
    al=np.asarray(al,float)
    d=cat.thx.astype(float)*(al[cat.a,0]-al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]-al[cat.b,1])
    s=cat.thx.astype(float)*(al[cat.a,0]+al[cat.b,0])+cat.thy.astype(float)*(al[cat.a,1]+al[cat.b,1])
    return d,s


def test_derivative_map_reduces_to_scalar_coefficient():
    """dalpha = g alpha must reproduce the scalar-g1 fit exactly (score, mean field, response, amplitudes)."""
    cat,t=synthetic(); g=.004; reg=np.array([0,1,2])
    scalar=amplitude(cat,with_junk(t),g,reg)
    mapped=[Template(u.alpha,u.name,u.kind,dalpha=g*u.alpha) for u in with_junk(t)]
    r=amplitude(cat,mapped,0.0,reg)          # the fallback coefficient is irrelevant when every template has a map
    assert np.allclose(r.q,scalar.q,rtol=1e-12) and np.allclose(r.mf,scalar.mf,rtol=1e-12)
    assert np.allclose(r.F,scalar.F,rtol=1e-12) and np.allclose(r.A,scalar.A,rtol=1e-6)   # dalpha is stored in float32
    assert np.allclose(r.jk_samples,scalar.jk_samples,rtol=1e-6)


def test_derivative_map_formulas_per_template():
    """Independent derivative maps per template: the documented contraction of the eleven accumulators."""
    cat,t=synthetic(); reg=np.array([0,1,2]); rng=np.random.default_rng(5)
    ts=with_junk(t); maps=[rng.normal(size=u.alpha.shape)*1e-3 for u in ts]
    tm=[Template(u.alpha,u.name,u.kind,dalpha=m) for u,m in zip(ts,maps)]
    r=amplitude(cat,tm,0.0,reg); x=cat.accum.sum(axis=2)
    D=[pair_scalar_pairs(cat,u.alpha) for u in tm]; Dp=[pair_scalar_pairs(cat,m) for m in maps]
    for i in range(3):
        d,_=D[i]; dp,sp=Dp[i]
        assert np.isclose(r.q[i],np.sum(d*x[:,0]+dp*x[:,1]+sp*x[:,2]),rtol=1e-12)
        assert np.isclose(r.mf[i],np.sum(d*x[:,8]+dp*x[:,9]+sp*x[:,10]),rtol=1e-12)
        for j in range(3):
            dj,_=D[j]; dpj,spj=Dp[j]
            F=np.sum(d*dj*x[:,3]+(d*dpj+dp*dj)*x[:,4]+dp*dpj*x[:,5]+.5*(d*spj+dj*sp)*x[:,6]+sp*spj*x[:,7])
            assert np.isclose(r.F[i,j],F,rtol=1e-12)
    assert np.allclose(r.F,r.F.T)


def test_joint_partials_match_single_fit():
    """joint_fit.partials_by_region (matrix products per region) equals amplitude._partials (per-pair arrays)."""
    from lyalenser.amplitude import _partials
    from lyalenser.joint_fit import partials_by_region
    cat,t=synthetic(); reg=np.array([0,1,0]); rng=np.random.default_rng(7)
    tm=[Template(u.alpha,u.name,u.kind,dalpha=rng.normal(size=u.alpha.shape)*1e-3) for u in with_junk(t)]
    tm[1]=Template(tm[1].alpha,tm[1].name,tm[1].kind)      # one template on the scalar fallback
    n1,r1,q1,F1,m1=_partials(cat,tm,.002,reg); n2,r2,q2,F2,m2=partials_by_region(cat,tm,.002,reg)
    assert n1==n2 and np.array_equal(r1,r2)
    assert np.allclose(q1,q2,rtol=1e-12) and np.allclose(F1,F2,rtol=1e-12) and np.allclose(m1,m2,rtol=1e-12)


def test_derivative_ratio_and_scalar_copies():
    from lyalenser.templates import derivative_ratio, with_scalar_derivative
    cat,t=synthetic(); ts=with_junk(t)
    tm=[Template(u.alpha,u.name,u.kind,dalpha=.5*u.alpha) for u in ts]
    assert np.isclose(derivative_ratio(tm),.5)
    copies=with_scalar_derivative(ts,.25)
    assert all(np.allclose(c.dalpha,.25*u.alpha) for c,u in zip(copies,ts)) and all(u.dalpha is None for u in ts)
