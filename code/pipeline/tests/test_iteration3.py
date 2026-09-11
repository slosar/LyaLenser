import numpy as np
from scipy.integrate import trapezoid
import pytest
from config import Config,SightlineSet
from mock import project_continuum,_interp3_chunked,rsd_displacement,sky_rays,ray_points,_deflection,sample_lognormal_quasars
from templates import flat_sky_band_filters,flat_sky_band_templates,Template,matched_shot_noise,matched_template_flat
from amplitude import amplitude
from validation_stats import absolute_statistics,slope_statistics,hotelling_shape
from pairs import PairCatalogue,find_pairs,accumulate,accumulate_slabs,pair_geometry
from xi_model import xi_from_data,xi_from_model,XiTable
from test_pairs import sample,table

def test_weighted_continuum_projection():
    rng=np.random.default_rng(100); x=np.arange(71)+3500.
    v=rng.normal(size=(5,71))+4+.1*x; w=rng.uniform(.1,4,v.shape)
    r=project_continuum(v,x,w)
    assert np.max(abs((r*w).sum(axis=1)))<1e-10
    assert np.max(abs((r*w*(x-x.mean())).sum(axis=1)))<1e-8
    assert np.allclose(project_continuum(r,x,w),r,atol=1e-12)

def test_periodic_interpolation_boundary():
    g=np.arange(4,dtype='f4').reshape(4,1,1)
    assert _interp3_chunked(g,np.array([[3.5,0.,0.],[-.5,0.,0.]]))==pytest.approx([1.5,1.5])

def test_fft_deflection_through_both_implementations():
    from templates import _flat_alpha_grid
    n=128; pix=2*np.pi/n; L=7.; x=np.arange(n)*pix
    phi0=2e-6; k=.5*L*L*phi0*np.cos(L*x)[:,None]*np.ones((1,n))
    target=-phi0*L*np.sin(L*x)[:,None]*np.ones((1,n))
    for a,b in (_flat_alpha_grid(k,pix),_deflection(k,pix,1.)):
        assert np.allclose(a,target,atol=1e-11); assert np.max(abs(b))<1e-12

def test_complete_band_unit_unfiltered_recovery():
    n=128; pix=2*np.pi/(n*5); ell,filters=flat_sky_band_filters((n,n),pix)
    assert np.max(abs(sum(filters.values())-1))<1e-15
    rng=np.random.default_rng(101); q=1800
    ra=180+np.rad2deg(rng.uniform(-n*pix/2,n*pix/2,q))/np.cos(np.deg2rad(30)); dec=30+np.rad2deg(rng.uniform(-n*pix/2,n*pix/2,q))
    x=np.arange(n)*pix; y=x.copy(); k=np.zeros((n,n))
    for L in (20,40,45,50,90,95,100,105,110,150,190,195,200,205,250,290,295,300,330):
        k+=np.cos(L*x[:,None]+.3)+.7*np.sin(L*y[None,:]+.8)
    basis,_=flat_sky_band_templates(k,ra,dec,pix)
    a=np.arange(q//2); b=a+q//2; tx=rng.normal(size=len(a)); ty=rng.normal(size=len(a)); norm=np.hypot(tx,ty); tx/=norm; ty/=norm
    acc=np.zeros((len(a),11,6)); acc[:,3,0]=1
    cat=PairCatalogue(a,b,tx,ty,np.ones(len(a)),acc,np.ones(len(a)))
    from amplitude import pair_scalars
    d,_,_=pair_scalars(cat,basis); signal=d[:3].sum(axis=0)+d[-1]
    cat.accum[:,0,0]=signal
    r=amplitude(cat,basis)
    assert r.A==pytest.approx([1,1,1,0,0,0,1],abs=1e-10)
    for t in basis[:3]: t.kind='truth'
    assert amplitude(cat,basis).A==pytest.approx([1,1,1,0,0,0,1],abs=1e-10)
    with pytest.raises(ValueError,match='missing required'): amplitude(cat,basis[:-2]+basis[-1:])
    broken=list(basis); broken[-1]=Template(basis[0].alpha,'junk','junk')
    with pytest.raises(ValueError,match='singular'): amplitude(cat,broken)

def test_spherical_rays_match_pair_metric_at_patch_positions():
    for ra in (171,180,189):
        for dec in (21,30,39):
            ras=np.array([ra,ra+.01]); decs=np.array([dec,dec+.012])
            _,_,theta=pair_geometry(ras,decs,np.array([0]),np.array([1]))
            for chi in (3500,4000,4500):
                points=sky_rays(ras,decs)*chi
                assert np.linalg.norm(points[0]-points[1])/(chi*theta[0])==pytest.approx(1,rel=1e-6)
                # All field consumers use the same sky->box function.
                grid=ray_points(ras,decs,chi,2,.5,(100,100,100),3000)
                assert np.allclose((grid[0]-grid[1])*[2,2,.5],points[0]-points[1])

def test_rsd_full_mode_spectrum_and_no_aliasing():
    n=32; dz=.5; z=np.arange(n)*dz; kz=2*np.pi*11/(n*dz)
    delta=np.broadcast_to(np.cos(kz*z),(16,16,n)).copy()
    v=rsd_displacement(delta,2,dz)
    truth=-.97/kz*np.sin(kz*z)
    assert np.allclose(v,truth,atol=1e-12)
    power=abs(np.fft.rfft(v[0,0]))**2
    assert power[11]/power.sum()==pytest.approx(1,abs=1e-12)

def test_integrated_poisson_intensity():
    g=np.full((4,4,4),.2,np.float32); rng=np.random.default_rng(102)
    q=sample_lognormal_quasars(g,2,.5,3500,100,1000,rng,angular_weight=np.full((4,4),.5))
    expected=100000*np.exp(3.5*.2)*.5
    assert q['integrated_intensity']==pytest.approx(expected,rel=1e-6)
    assert abs(len(q['x'])-expected)<5*np.sqrt(expected)

def test_shot_noise_uniform_reduction():
    from cosmo import chi,z_of_chi
    from cross_spectrum import kernel,Z_CMB
    edges=np.linspace(3500,4500,41); n2d=1e5; D=1000.; b=3.5; ratio=.05
    shot=matched_shot_noise(edges,np.full(40,n2d/D),lambda z:np.full_like(z,b),ratio)
    c=np.linspace(3500,4500,10001)
    expected=(1+ratio)*D*trapezoid(kernel(c,float(chi(Z_CMB)))**2,c)/(b*b*n2d)
    assert shot==pytest.approx(expected,rel=1e-8)

def test_slab_midpoint_and_invalid_pixels():
    from cosmo import z_of_chi
    sl=sample(4,12); boundary=3904.
    sl.slab[:]=(sl.chi>=boundary).astype('i1')
    cfg=Config(chi_ref=3900,slabs=((float(z_of_chi(3890)),float(z_of_chi(boundary))),(float(z_of_chi(boundary)),float(z_of_chi(3920)))))
    pairs=find_pairs(sl,30/sl.chi.min()); cats=accumulate_slabs(sl,pairs,table(),cfg)
    whole=accumulate(sl,pairs,table(),cfg)
    assert sum(c.npair.sum() for c in cats.values())==whole.npair.sum()
    assert all(c.npair.sum()>0 for c in cats.values())
    sl.slab[:]=-1
    assert len(accumulate(sl,pairs,table(),cfg).a)==0
    assert xi_from_data(sl,cfg).meta['accepted_weight']==0

def test_absolute_statistics_do_not_cancel_contamination():
    x=np.array([0,1,5,10]); y=7+np.arange(20)[:,None]*.01+x
    assert slope_statistics(y)['mean']==pytest.approx(1)
    stats=absolute_statistics(y[:,1],1)
    assert stats['residual']>7 and stats['bound95']>7
    rng=np.random.default_rng(103); assert 0<=hotelling_shape(rng.normal(size=(20,6)))['p_value']<=1

def test_amplitude_interpolation_convergence():
    # Actual compression and matrix fit, evaluated off both grids.
    rng=np.random.default_rng(104); sl=sample(20,30); cfg=Config(chi_ref=3900)
    from templates import curl
    signals=[Template(rng.normal(size=(20,2)),str(lo),'signal',Lmin=lo,Lmax=hi) for lo,hi in ((40,100),(100,200),(200,300))]
    basis=signals+[Template(curl(t.alpha),t.name+'_curl','curl',Lmin=t.Lmin,Lmax=t.Lmax) for t in signals]+[Template(rng.normal(size=(20,2)),'junk','junk')]
    values=[]
    for step in (.25,.125):
        grid=np.arange(0,40+step/2,step); rp=grid[:,None]; rz=grid[None,:]
        xi=np.exp(-(rp*rp+rz*rz)/200)
        tab=XiTable(grid,grid,xi,-rp/100*xi)
        cat=accumulate(sl,find_pairs(sl,30/sl.chi.min()),tab,cfg)
        values.append(amplitude(cat,basis).A)
    assert np.linalg.norm(values[1]-values[0])/np.linalg.norm(values[1])<.005

def test_analytic_derivative_nk_convergence():
    from forest_power import ForestPower
    cfg=Config(xi_step=2,xi_max=30)
    pf=ForestPower(model='kaiser')
    a=xi_from_model(pf,cfg,nk=6400); b=xi_from_model(pf,cfg,nk=12800)
    take=a.r_perp>=10
    assert np.linalg.norm(a.xi_rp[take]-b.xi_rp[take])/np.linalg.norm(b.xi_rp[take])<.01

def test_matched_flat_and_spherical_end_to_end():
    import healpy as hp
    from templates import matched_template
    from cosmo import z_of_chi
    rng=np.random.default_rng(104); nside=32
    ra,dec=hp.pix2ang(nside,np.arange(hp.nside2npix(nside)),lonlat=True)
    footprint=(abs((ra-180)*np.cos(np.deg2rad(30)))<10)&(abs(dec-30)<10)
    area=footprint.sum()*hp.nside2pixarea(nside); side=np.sqrt(area)
    n=10000; r=180+np.rad2deg(rng.uniform(-side*.4,side*.4,n))/np.cos(np.deg2rad(30)); d=30+np.rad2deg(rng.uniform(-side*.4,side*.4,n))
    c=rng.uniform(3500,4500,n); q={'ra':r,'dec':d,'z':z_of_chi(c)}
    bias=lambda z:np.full_like(z,3.5)
    flat,mask,meta=matched_template_flat(q,q,bias,(32,32),side/32,nmin_rand=0)
    alm,mask_s,meta_s=matched_template(q,q,bias,Config(nside_alpha=nside,lmax_alpha=32),nmin_rand=0,footprint_mask=footprint)
    assert np.max(abs(flat))==0 and np.max(abs(alm))==0
    assert meta['shot_s']==pytest.approx(meta_s['shot_s'],rel=1e-10)
    assert meta['completeness'].mean()==pytest.approx(1)
    assert meta_s['completeness'][footprint].mean()==pytest.approx(1)
    # A spatial density modulation must survive radial matched weighting.
    keep=rng.random(n)<(.5+.4*np.cos(np.deg2rad((r-180)*np.cos(np.deg2rad(30)))*20))
    data={k:v[keep] for k,v in q.items()}
    fm,fmask,_=matched_template_flat(data,q,bias,(32,32),side/32,nmin_rand=0)
    _,smask,smeta=matched_template(data,q,bias,Config(nside_alpha=nside,lmax_alpha=32),nmin_rand=0,footprint_mask=footprint)
    x=(np.arange(32)+.5-16)*side/32
    fcos=np.broadcast_to(np.cos(20*x)[:,None],fm.shape)
    scos=np.cos(20*np.deg2rad((ra-180)*np.cos(np.deg2rad(30))))
    def fit_mode(field,mode,mask):
        u=mode[mask]-mode[mask].mean()
        return np.sum(field[mask]*u)/np.sum(u*u)
    fa=fit_mode(fm,fcos,fmask); sa=fit_mode(smeta["kappa_map"],scos,smask)
    assert fa>0 and sa==pytest.approx(fa,rel=.1)

def test_pixel_prediction_matches_dense_projector():
    from response import prediction_catalogue
    from amplitude import pair_scalars
    rng=np.random.default_rng(100); sl=sample(6,15); cfg=Config(chi_ref=3900)
    tab=table(); cat=accumulate(sl,find_pairs(sl,30/sl.chi.min()),tab,cfg)
    mod=rng.normal(scale=.2,size=len(sl.chi)); pred=prediction_catalogue(sl,cat,mod,tab,tab,cfg)
    for ip,(a,b) in enumerate(zip(cat.a,cat.b)):
        ia=slice(sl.pix_start[a],sl.pix_start[a+1]); ib=slice(sl.pix_start[b],sl.pix_start[b+1])
        ca=sl.chi[ia].astype(float); cb=sl.chi[ib].astype(float); wa=sl.w[ia].astype(float); wb=sl.w[ib].astype(float)
        def P(c,w):
            u=np.column_stack((np.ones(len(c)),c-c.mean()))
            return np.eye(len(c))-u@np.linalg.solve(u.T@(w[:,None]*u),u.T*w)
        rp=(ca[:,None]+cb)*.5*cat.theta[ip]; rz=abs(ca[:,None]-cb)
        C=np.array([[tab.interp(x,y)[0] for x,y in zip(xx,yy)] for xx,yy in zip(rp,rz)])
        B=C*((1+mod[ia,None])*(1+mod[ib][None,:])-1)
        expected=P(ca,wa)@B@P(cb,wb).T
        G=np.array([[tab.interp(x,y)[1] for x,y in zip(xx,yy)] for xx,yy in zip(rp,rz)])*(ca[:,None]+cb)*.5
        assert pred.accum[ip,0].sum()==pytest.approx(np.sum(wa[:,None]*wb*expected*G),rel=1e-5,abs=1e-7)

def test_cross_boundary_pixel_pair_is_assigned_by_midpoint():
    from cosmo import chi
    zedge=2.45; c=float(chi(zedge))
    sl=SightlineSet([0,1],[180,180.01],[30,30],[3,3],[0,1,2],[c-1,c+2],[1,1],[1,1],[0,1])
    cfg=Config(slabs=((2.1,zedge),(zedge,3.0)))
    cats=accumulate_slabs(sl,find_pairs(sl,30/sl.chi.min()),table(),cfg)
    assert cats['slab0'].npair.sum()==0
    assert cats['slab1'].npair.sum()==1

def test_disjoint_template_sightline_selection():
    from mock import generate_mock
    m=generate_mock(Config(scale=.03),seed=104,scale=.03,response=False,disjoint_selection=True,pixel_noise_power=0)
    assert len(m.sightlines.qid)>0 and len(m.quasars['qid'])>0
    assert len(np.intersect1d(m.sightlines.qid,m.quasars['qid']))==0

def test_cross_field_projection_at_multiple_patch_positions():
    from mock import project_lightcone_fields,ray_points,_interp3_chunked
    # A known periodic plane wave, spanning a 20-degree field. Independent
    # direct trilinear samples verify each projection and the QSO intensity.
    n=32; nz=48; cref=4000.; dx=cref*np.deg2rad(20)/n; dz=20.
    chis=3500+np.arange(nz)*dz; x=np.arange(n)*2*np.pi/n
    dm=np.broadcast_to(np.cos(3*x)[:,None,None],(n,n,nz)).astype('f4').copy()
    long=.2*dm; wc=np.linspace(.001,.002,nz); wl=np.maximum(0.,.003-(chis-3500)*1e-5); trap=np.full(nz,dz); trap[[0,-1]]*=.5
    intensity,kc,kl,dl=project_lightcone_fields(dm,long,chis,dx,dz,cref,wc,wl,trap,(3600,4200))
    for i,j in ((4,7),(16,16),(27,24)):
        ra=180+np.rad2deg((i-n/2)*dx/cref)/np.cos(np.deg2rad(30)); dec=30+np.rad2deg((j-n/2)*dx/cref)
        pts=ray_points(np.full(nz,ra),np.full(nz,dec),chis,dx,dz,dm.shape,chis[0])
        sampled=_interp3_chunked(dm,pts)
        assert np.allclose(intensity[i,j],sampled,atol=1e-7)
        assert kc[i,j]==pytest.approx(np.dot(wc*trap,sampled),rel=2e-6,abs=1e-7)
        assert kl[i,j]==pytest.approx(np.dot(wl*trap,sampled),rel=2e-6,abs=1e-7)
        mask=(chis>=3600)&(chis<=4200)
        assert dl[i,j]==pytest.approx(.2*sampled[mask].mean(),abs=1e-7)

def test_hotelling_uses_full_covariance_and_f_calibration():
    from scipy.stats import f
    rng=np.random.default_rng(104); transform=np.tril(np.full((6,6),.7))+np.eye(6)
    values=rng.normal(size=(40,6))@transform.T+np.linspace(0,.7,6)
    differences=values[:,1:]-values[:,0,None]
    mu=differences.mean(axis=0); covariance=np.cov(differences.T)
    expected=40*mu@np.linalg.inv(covariance)@mu
    result=hotelling_shape(values)
    assert result['T2']==pytest.approx(expected,rel=1e-12)
    assert result['p_value']==pytest.approx(f.sf(35*expected/(5*39),5,35),rel=1e-12)
    diagonal=40*np.sum(mu**2/np.diag(covariance))
    assert abs(diagonal-expected)>.1

def test_exact_trilinear_covariance_window():
    from grid_covariance import sample_covariance,b3
    # White independent grid pixels: phase-averaged variance is (2/3)^3.
    cube=np.zeros((13,13,13)); cube[6,6,6]=1
    grid=np.array([0.,.25,.5,1.]); x,g=sample_covariance(cube,grid,1.,1.,128)
    assert x[0,0]==pytest.approx((2/3)**3,abs=1e-14)
    assert x[0,2]==pytest.approx((2/3)**2*b3(.5)[0],abs=1e-14)
    assert abs(g[0,0])<1e-14
    x2,g2=sample_covariance(cube,grid,1.,1.,256)
    assert np.allclose(x,x2,atol=1e-7)
    assert np.allclose(g,g2,atol=1e-6)
