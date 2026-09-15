"""Low-redshift tracers for the forest-lensing cross-correlation (Stage A iteration 7).

Why: the CMB-lensing cross-correlation measures lensing plus the intrinsic forest-forest-density three-point
function of the forest's own redshift range (NOTES.md iteration 6; CMB_FUTURE_WORK.md). A density tracer at
z < 1.8 shares no modes with the forest at z = 2.1-3, so its correlation with the forest pair products is lensing
only (plus magnification of the sightline quasars, handled by the mean field, and shared sky systematics).

What this module provides
- A tracer table (name, redshift slice, bias, surface density; DR1-like defaults for the mock; the data path
  measures n(z) from the catalogues and b from the angular auto-spectrum).
- Limber spectra per redshift slice: the 3 x 3 covariance of (slice-averaged matter density, the slice's
  contribution to kappa_lya, its contribution to kappa_CMB), and the spectra of the foreground not covered by
  slices plus the background of kappa_CMB. Slices are independent in the Limber approximation.
- Mock realisation: correlated Gaussian maps per slice, 2-D lognormal Poisson tracers per (tracer, slice), the
  rest of the foreground/background as before; kappa_lya_rest = sum of slice contributions + rest.
- Templates per seed: kernel-weighted tracer maps (forest-source lensing kernel W(chi; chi_ref), per-object
  1/(b nbar)), the bias b of every tracer from its own masked angular auto-spectrum (shot noise subtracted, mask
  coupling and pixel window applied to the theory), the Wiener combination of the tracers of a slice into one
  estimate of the slice's kappa_lya contribution, and the sum over slices as the combined estimate (slices are
  independent, so the sum of the Wiener-filtered slice maps is the conditional expectation of kappa_lya given
  all tracers); band templates (3 science bands + 3 curl partners + junk) for every slice and for the sum.
- Amplitudes: one fit per slice and one for the combined template with the same jackknife regions, the joint
  jackknife covariance of the slice amplitudes, and their optimal combination as a cross-check of the
  combined-template amplitude. Sightline pixels carry the inverse-variance weights w = 1/(sigma_N^2 + sigma_F^2)
  (pixel diagonal of C^-1; correlations along the skewer are not included in the weights).
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from functools import lru_cache
from pathlib import Path
import sys, json
import numpy as np
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
if str(CODE) not in sys.path: sys.path.insert(0,str(CODE))
from cosmo import chi as chi_of_z, z_of_chi
from cross_spectrum import kernel, Z_CMB
from recon_noise import DEG2

SCIENCE_L=(40.,300.)
ANNULUS=40.
KMAX_BIAS=0.2        # h/Mpc: the bias is fitted on large scales only, ell <= KMAX_BIAS * chi(z_mid) of the slice (user, 2026-09-15)
LMIN_BIAS=40.


@dataclass(frozen=True)
class Tracer:
    name: str
    zmin: float
    zmax: float
    bias: float          # generator bias (mock); the pipeline measures b from the auto-spectrum
    nbar_deg2: float     # surface density in the slice (mock)
    smoothing_mpc: float = 3.0   # transverse Gaussian smoothing (comoving, at the slice's mid distance) of the projected
                                 # density before the lognormal transform: without it the 0.7 Mpc/h generator cells give
                                 # b^2 var(delta_2D) ~ 2.6 in the lowest slice and the tracer decorrelates from the density.
                                 # The filter G(l) = exp(-l^2 sigma^2 / 2), sigma = smoothing / chi_mid, is a known part of
                                 # the mock tracer model and enters the template theory (data path: no smoothing).

    @property
    def chi_mid(self): return float(chi_of_z(.5*(self.zmin+self.zmax)))
    def filter(self,ell):
        sig=self.smoothing_mpc/self.chi_mid
        return np.exp(-.5*(np.asarray(ell,float)*sig)**2)

    @property
    def label(self): return f'{self.name}_{self.zmin:g}_{self.zmax:g}'
    @property
    def chi_range(self): return (float(chi_of_z(self.zmin)),float(chi_of_z(self.zmax)))


# DR1-like mock inputs. Slices end at z = 1.75, below the box front (chi(2.1) - 300 Mpc/h, z = 1.805): no tracer
# shares a density mode with the forest box. b_QSO(z) from the DESI fit 0.278((1+z)^2-6.565)+2.393.
TRACERS=(Tracer('LRG',.4,.6,1.9,200.),Tracer('LRG',.6,.8,2.1,250.),Tracer('LRG',.8,1.1,2.3,150.),
         Tracer('ELG',.8,1.1,1.3,350.),Tracer('ELG',1.1,1.6,1.4,500.),
         Tracer('QSO',.8,1.1,1.7,40.),Tracer('QSO',1.1,1.6,2.1,60.),Tracer('QSO',1.6,1.75,2.5,25.))


def slices_of(tracers=TRACERS):
    """Sorted redshift slices and the tracers in each; every tracer must cover exactly one slice."""
    bounds=sorted({(t.zmin,t.zmax) for t in tracers})
    for (a,b),(c,d) in zip(bounds[:-1],bounds[1:]):
        if c<b: raise ValueError(f'overlapping slices {(a,b)} and {(c,d)}')
    return [{'zmin':a,'zmax':b,'tracers':[t for t in tracers if (t.zmin,t.zmax)==(a,b)]} for a,b in bounds]


def _ell_grid(lmax):
    return np.unique(np.r_[2.,np.linspace(2,max(2.,lmax),640)])


@lru_cache(maxsize=32)
def slice_spectra(zmin,zmax,chi_source,lmax):
    """Limber 3 x 3 covariance per L of (slice-mean density d, kappa_lya part l, kappa_CMB part c) over [zmin, zmax]."""
    from three_tracer import limber
    c1,c2=float(chi_of_z(zmin)),float(chi_of_z(zmax)); ccmb=float(chi_of_z(Z_CMB))
    wd=lambda c: np.where((c>=c1)&(c<=c2),1./(c2-c1),0.)
    wl=lambda c: kernel(c,chi_source); wc=lambda c: kernel(c,ccmb)
    L=_ell_grid(lmax); W=(wd,wl,wc); C=np.zeros((len(L),3,3))
    for i in range(3):
        for j in range(i,3):
            C[:,i,j]=C[:,j,i]=limber(L,W[i],W[j],c1,c2,nchi=400,to_recombination=False)
    return L,C


@lru_cache(maxsize=32)
def rest_spectra(front_ranges,cbox_max,chi_source,lmax):
    """(ll, lc, cc) of the foreground ranges not covered by slices, plus cc behind the box."""
    from three_tracer import limber
    ccmb=float(chi_of_z(Z_CMB)); wl=lambda c: kernel(c,chi_source); wc=lambda c: kernel(c,ccmb)
    L=_ell_grid(lmax); ll=np.zeros(len(L)); lc=np.zeros(len(L)); cc=np.zeros(len(L))
    for a,b in front_ranges:
        if b-a<1.: continue
        ll+=limber(L,wl,wl,a,b,nchi=400,to_recombination=False)
        lc+=limber(L,wl,wc,a,b,nchi=400,to_recombination=False)
        cc+=limber(L,wc,wc,a,b,nchi=400,to_recombination=False)
    if cbox_max<ccmb: cc+=limber(L,wc,wc,cbox_max,ccmb,nchi=500,to_recombination=True)
    return L,ll,lc,cc


def front_ranges(slices,cbox_min):
    """Foreground chi ranges between 1 Mpc/h and the box front that no slice covers."""
    out=[]; start=1.
    for s in slices:
        c1,c2=float(chi_of_z(s['zmin'])),float(chi_of_z(s['zmax']))
        if c2>cbox_min: raise ValueError(f"slice {s['zmin']}-{s['zmax']} reaches into the forest box (front at chi={cbox_min:.0f})")
        out.append((start,c1)); start=c2
    out.append((start,float(cbox_min)))
    return tuple((a,b) for a,b in out if b-a>1.)


def gaussian_maps(nx,ny,side_angle,L,cov,rng):
    """n correlated Gaussian flat-sky maps from a per-multipole covariance cov[L, n, n] (Cholesky per mode)."""
    lx=2*np.pi*np.fft.fftfreq(nx,side_angle/nx); ly=2*np.pi*np.fft.rfftfreq(ny,side_angle/ny)
    ell=np.hypot(lx[:,None],ly[None,:]); n=cov.shape[1]; pixarea=(side_angle/nx)*(side_angle/ny)
    C=np.empty(ell.shape+(n,n))
    for i in range(n):
        for j in range(n):
            C[...,i,j]=np.interp(ell,L,cov[:,i,j],left=0.,right=cov[-1,i,j])
    C[...,range(n),range(n)]+=1e-12*np.max(np.abs(cov))
    root=np.linalg.cholesky(C)
    z=np.stack([np.fft.rfft2(rng.normal(size=(nx,ny)).astype(np.float32)) for _ in range(n)],axis=-1)
    f=np.einsum('...ij,...j->...i',root,z)/np.sqrt(pixarea); f[0,0,:]=0
    return [np.fft.irfft2(f[...,i],s=(nx,ny)).astype(np.float32) for i in range(n)]


def sample_lognormal_2d(delta,bias,nbar_deg2,side_angle,rng,angular_weight=None):
    """Poisson sample of the 2-D lognormal intensity nbar exp(b delta - b^2 var/2) x completeness.
    Returns angular offsets (rad) from the patch centre; the mean count is nbar x area x <completeness>."""
    nx,ny=delta.shape; var=float(np.var(delta)); corr=.5*bias*bias*var
    w=np.exp(np.clip(bias*np.asarray(delta,np.float64)-corr,-30,30))
    if angular_weight is not None: w=w*np.asarray(angular_weight,np.float64)
    area_deg2=np.rad2deg(side_angle)**2; lam=nbar_deg2*area_deg2/(nx*ny)*w
    counts=rng.poisson(lam); ix,iy=np.nonzero(counts); rep=counts[ix,iy]
    ix=np.repeat(ix,rep); iy=np.repeat(iy,rep); n=len(ix); pix=side_angle/nx
    # Node-centred periodic cells, as the 3-D quasar sampler: cell i covers [i-1/2, i+1/2) around the field node.
    x=((ix+rng.random(n)-.5)%nx-nx/2)*pix; y=((iy+rng.random(n)-.5)%ny-ny/2)*pix
    return {'x':x,'y':y,'ix':ix,'iy':iy,'density_variance':var,'lognormal_correction':corr,'expected':float(lam.sum())}


def angles_to_radec(x,y,center=(180.,30.)):
    return center[0]+np.rad2deg(x)/np.cos(np.deg2rad(center[1])),center[1]+np.rad2deg(y)


def lowz_realisation(nx,ny,side_angle,cbox,chi_source,streams,completeness=None,tracers=TRACERS,random_factor=20):
    """Foreground of the mock: per-slice (density, kappa_lya part, kappa_CMB part) Gaussian maps, lognormal Poisson
    tracers, the uncovered rest of the foreground and the background of kappa_CMB.
    Returns rest_lya, rest_cmb (to be added to the box contributions), maps, catalogue, randoms, spectra, attrs."""
    slices=slices_of(tracers)
    lx=2*np.pi*np.fft.fftfreq(nx,side_angle/nx); ly=2*np.pi*np.fft.rfftfreq(ny,side_angle/ny)
    lmax=int(np.ceil(np.hypot(lx[:,None],ly[None,:]).max()))
    rng_f=streams['lowz_fields']; rng_s=streams['lowz_sampling']; rng_r=streams['lowz_randoms']
    maps={}; rest_lya=np.zeros((nx,ny),np.float32); rest_cmb=np.zeros((nx,ny),np.float32); spectra={'slices':[]}
    cat={k:[] for k in ('ra','dec','z','chi','tracer','slice')}; rnd={k:[] for k in ('ra','dec','z','chi','tracer','slice')}
    for i,s in enumerate(slices):
        L,C=slice_spectra(s['zmin'],s['zmax'],float(chi_source),lmax)
        d,kl,kc=gaussian_maps(nx,ny,side_angle,L,C,rng_f)
        maps[f'lowz_delta_{i}']=d; maps[f'lowz_kappa_lya_{i}']=kl; maps[f'lowz_kappa_cmb_{i}']=kc
        rest_lya+=kl; rest_cmb+=kc
        c1,c2=float(chi_of_z(s['zmin'])),float(chi_of_z(s['zmax']))
        spectra['slices'].append({'zmin':s['zmin'],'zmax':s['zmax'],'chi':[c1,c2],'L':L,'C':C})
        for t in s['tracers']:
            ti=tracers.index(t)
            if t.smoothing_mpc>0:
                from scipy.ndimage import gaussian_filter
                ds=gaussian_filter(d,sigma=t.smoothing_mpc/(t.chi_mid*side_angle/nx),mode='wrap')
            else: ds=d
            g=sample_lognormal_2d(ds,t.bias,t.nbar_deg2,side_angle,rng_s,completeness)
            ra,dec=angles_to_radec(g['x'],g['y']); chi=rng_s.uniform(c1,c2,len(ra))
            cat['ra'].append(ra); cat['dec'].append(dec); cat['chi'].append(chi); cat['z'].append(z_of_chi(chi))
            cat['tracer'].append(np.full(len(ra),ti,np.int16)); cat['slice'].append(np.full(len(ra),i,np.int16))
            nr=random_factor*len(ra)
            if completeness is None:
                rx=rng_r.uniform(-side_angle/2,side_angle/2,nr); ry=rng_r.uniform(-side_angle/2,side_angle/2,nr)
            else:
                from mock import _sample_angular_selection
                rx,ry=_sample_angular_selection(nr,rng_r,side_angle,np.asarray(completeness)/float(np.max(completeness)))
            rra,rdec=angles_to_radec(rx,ry); rchi=rng_r.uniform(c1,c2,nr)
            rnd['ra'].append(rra); rnd['dec'].append(rdec); rnd['chi'].append(rchi); rnd['z'].append(z_of_chi(rchi))
            rnd['tracer'].append(np.full(nr,ti,np.int16)); rnd['slice'].append(np.full(nr,i,np.int16))
    L,ll,lc,cc=rest_spectra(front_ranges(slices,cbox[0]),float(cbox[1]),float(chi_source),lmax)
    cov=np.zeros((len(L),2,2)); cov[:,0,0]=ll; cov[:,0,1]=cov[:,1,0]=lc; cov[:,1,1]=cc
    rl,rc=gaussian_maps(nx,ny,side_angle,L,cov,rng_f)
    maps['lowz_rest_lya']=rl; maps['lowz_rest_cmb']=rc; rest_lya+=rl; rest_cmb+=rc
    spectra['rest']={'L':L,'klkl_rest':ll,'klkc_rest':lc,'kckc_rest':cc,'front_ranges':front_ranges(slices,cbox[0])}
    # Totals of the foreground/background, for the spectrum check (same keys as the CMB-only generator).
    tot={'L':L,'klkl_rest':ll.copy(),'klkc_rest':lc.copy(),'kckc_rest':cc.copy()}
    for s in spectra['slices']:
        tot['klkl_rest']+=np.interp(L,s['L'],s['C'][:,1,1]); tot['klkc_rest']+=np.interp(L,s['L'],s['C'][:,1,2]); tot['kckc_rest']+=np.interp(L,s['L'],s['C'][:,2,2])
    cat={k:np.concatenate(v) if v else np.zeros(0) for k,v in cat.items()}
    rnd={k:np.concatenate(v) if v else np.zeros(0) for k,v in rnd.items()}
    attrs={'tracers':[asdict(t) for t in tracers],'slices':[{'zmin':s['zmin'],'zmax':s['zmax'],'tracers':[tracers.index(t) for t in s['tracers']]} for s in slices],
           'random_factor':random_factor,'lognormal':'2-D projected slice density, Poisson sampled'}
    return rest_lya,rest_cmb,maps,cat,rnd,spectra,tot,attrs


# ----------------------------------------------------------------------------------------------------------------
# Templates from a tracer catalogue (mock: flat patch; the data path builds HEALPix maps with the same weights).

def _annuli(pix,n):
    lx=2*np.pi*np.fft.fftfreq(n,pix); ell=np.hypot(lx[:,None],lx[None,:])
    edges=np.arange(0,np.sqrt(2)*np.pi/pix+2*ANNULUS,ANNULUS); L=.5*(edges[1:]+edges[:-1])
    return ell,edges,L


def mask_coupled(theory_L,theory,ell,edges,mask,n,window=1.):
    """Binned pseudo-spectrum expected from a theory spectrum under the common mask (and pixel window)."""
    grid=np.interp(ell,theory_L,theory,left=0,right=0)*window
    maskpower=abs(np.fft.fft2(mask))**2
    convolved=np.fft.ifft2(np.fft.fft2(maskpower)*np.fft.fft2(grid)).real/n**4/max(np.mean(np.asarray(mask)**2),1e-30)
    count=np.histogram(ell,edges)[0]
    return np.divide(np.histogram(ell,edges,weights=convolved)[0],count,out=np.zeros(len(edges)-1),where=count>0)


def bias_band(zmin,zmax,kmax=KMAX_BIAS,lmin=LMIN_BIAS):
    """Multipole range of the bias fit: lmin <= ell <= kmax chi(z_mid) (Limber k = (ell + 1/2)/chi)."""
    cmid=float(chi_of_z(.5*(zmin+zmax)))
    return (float(lmin),float(kmax*cmid-.5))


def fit_bias(measured,shot,theory_binned,L,nmodes,band=SCIENCE_L):
    """b^2 = weighted least squares of (C_hat - shot) on the unit-bias theory over ``band`` (large scales only,
    bias_band); error from the Gaussian mode count. Returns b, sigma_b, b^2, chi2 of the fit."""
    use=(L>=band[0])&(L<=band[1])&(theory_binned>0)&(nmodes>0)
    # nmodes counts the full Fourier plane (+-l are conjugate): nmodes/2 independent modes per annulus.
    y=measured[use]-shot; T=theory_binned[use]; var=2*np.maximum(measured[use],1e-30)**2/(nmodes[use]/2)
    w=1/var; b2=float(np.sum(w*y*T)/np.sum(w*T*T)); vb2=1/float(np.sum(w*T*T))
    b=float(np.sqrt(max(b2,1e-12))); sb=float(np.sqrt(vb2)/(2*b))
    chi2=float(np.sum(w*(y-b2*T)**2))
    return {'b':b,'sigma_b':sb,'b2':b2,'sigma_b2':float(np.sqrt(vb2)),'chi2':chi2,'dof':int(use.sum()-1),'band':[float(band[0]),float(band[1])],
            'annuli_used':int(use.sum())}


def combine_maps(maps,weights,edges,pix):
    """Sum_k w_k(L) t_k in Fourier space with per-annulus weights (rows of `weights`: annuli; columns: maps)."""
    n=maps[0].shape[0]; lx=2*np.pi*np.fft.fftfreq(n,pix); ly=2*np.pi*np.fft.rfftfreq(n,pix)
    ell=np.hypot(lx[:,None],ly[None,:]); idx=np.clip(np.searchsorted(edges,ell,side='right')-1,0,len(edges)-2)
    out=np.zeros(ell.shape,complex)
    for k,m in enumerate(maps): out+=np.fft.rfft2(m)*weights[idx,k]
    return np.fft.irfft2(out,s=(n,n)).astype(np.float32)


def tracer_maps(mock,cfg,n=128,fixed_bias=None):
    """Per-tracer kernel-weighted maps, masks, shot noise, auto-spectra and fitted bias; per-slice Wiener combination;
    the combined map. Everything on the generator's angular patch side."""
    from templates import matched_template_flat, map_spectrum
    from mock import patch_side_rad
    side=patch_side_rad(mock); pix=side/n; shape=(n,n); cref=float(mock.sightlines.attrs['chi_ref'])
    ell,edges,L=_annuli(pix,n); nmodes=np.histogram(ell,edges)[0].astype(float)
    lx=2*np.pi*np.fft.fftfreq(n,pix); W1=np.sinc(lx[:,None]*pix/(2*np.pi))*np.sinc(lx[None,:]*pix/(2*np.pi)); window=W1**2   # NGP pixel window: W1 in a cross with the continuous field, W1^2 in the auto
    attrs=json.loads(mock.attrs['lowz']) if isinstance(mock.attrs.get('lowz'),str) else mock.attrs['lowz']
    tracers=[Tracer(**t) for t in attrs['tracers']]; slices=attrs['slices']
    cat=mock.lowz_catalogue; rnd=mock.lowz_randoms
    lmax=int(np.ceil(ell.max()))
    out={'pixel_size':pix,'side_rad':side,'edges':edges,'L':L,'nmodes':nmodes,'tracers':{},'slices':[],'window':window}
    combined=np.zeros(shape,np.float32); common_mask=np.ones(shape,np.float32)
    mock2d=str(attrs.get('lognormal','')).startswith('2-D')
    for si,s in enumerate(slices):
        Lth,C=slice_spectra(s['zmin'],s['zmax'],cref,lmax)
        # Signal model of a unit-bias kernel-weighted tracer map: for the mock's 2-D slice the map samples the single
        # projected density, so <t t> = (int W)^2 C_dd and <t kappa> = (int W) C_dl exactly; for real tracers (density
        # evolving along the slice) both are C_ll. The two differ by < 1.1 % here (NOTES iteration 7).
        if mock2d:
            c1,c2=float(chi_of_z(s['zmin'])),float(chi_of_z(s['zmax'])); cc=np.linspace(c1,c2,2000)
            Wint=float(np.trapz(kernel(cc,cref),cc)); S_auto=Wint*Wint*C[:,0,0]; S_cross=Wint*C[:,0,1]
        else: S_auto=C[:,1,1]; S_cross=C[:,1,1]
        S_ll=S_auto
        maps=[]; masks=[]; shots=[]; names=[]; filters=[]
        for ti in s['tracers']:
            t=tracers[ti]; sel=cat['tracer']==ti; rsel=rnd['tracer']==ti
            q={'ra':cat['ra'][sel],'dec':cat['dec'][sel],'z':cat['z'][sel]}; r={'ra':rnd['ra'][rsel],'dec':rnd['dec'][rsel],'z':rnd['z'][rsel]}
            m1,mask,meta=matched_template_flat(q,r,lambda z: np.ones_like(np.asarray(z,float)),shape,pix,
                                              radial_range=t.chi_range,source_chi=cref,radial_bins=1)
            # Unit-bias map: b G kappa_slice + noise (G = the mock tracer's smoothing filter, 1 for data). Its
            # auto-spectrum gives b^2 against the unit-bias theory C_ll G^2 W1^2 (mask-coupled).
            G=t.filter(ell); Ta=mask_coupled(Lth,S_ll,ell,edges,mask,n,window*G*G); auto=map_spectrum(m1,pix,edges,mask)
            fit=fit_bias(auto,meta['shot_s'],Ta,L,nmodes,band=bias_band(s['zmin'],s['zmax']))
            # A fit is usable when b^2 is positive and resolved (b^2 > 3 sigma); otherwise the generator's bias is
            # used and flagged (smoke patches have too few modes; at scale 1 no fallback is expected).
            usable=np.isfinite(fit['b2']) and fit['b2']>3*fit['sigma_b2']
            fit['usable']=bool(usable)
            if fixed_bias is not None: b=float(fixed_bias[t.label]); source='fixed'
            elif usable: b=fit['b']; source='auto-spectrum'
            else: b=float(t.bias); source='fallback_true'
            fit['bias_source']=source
            maps.append((m1/b).astype(np.float32)); masks.append(mask); shots.append(meta['shot_s']/b**2); names.append(t.label); filters.append(G)
            out['tracers'][t.label]={'slice':si,'bias_true':t.bias,'bias_fit':fit,'bias_used':b,'bias_source':source,'n_objects':int(sel.sum()),
                                    'shot_s_unit_bias':meta['shot_s'],'auto_spectrum':auto,'theory_unit_bias':Ta,'fsky':meta['fsky'],'smoothing_mpc':t.smoothing_mpc}
        mask=np.prod(masks,axis=0).astype(np.float32); common_mask*=mask
        # Wiener combination within the slice with the MODEL covariance: map k = G_k kappa_slice W1 + shot_k, so
        # C_kl = C_ll G_k G_l W1^2 (mask-coupled) + shot_k delta_kl and the signal vector <t_k kappa> = C_ll G_k W1.
        # w = C^-1 s: inverse-shot-noise weighting of the tracers, Wiener-suppressed where shot noise dominates,
        # and the known smoothing / pixel window deconvolved. Measured spectra are kept as a model diagnostic.
        k=len(maps); Ta=np.zeros((len(L),k,k)); Tc=np.zeros((len(L),k))
        for a in range(k):
            Tc[:,a]=mask_coupled(Lth,S_cross,ell,edges,mask,n,W1*filters[a])
            for c in range(a,k):
                Ta[:,a,c]=Ta[:,c,a]=mask_coupled(Lth,S_ll,ell,edges,mask,n,window*filters[a]*filters[c])
        Chat=np.zeros((len(L),k,k))
        for a in range(k):
            for c in range(a,k):
                Chat[:,a,c]=Chat[:,c,a]=map_spectrum(maps[a],pix,edges,mask,map_b=maps[c])
        weights=np.zeros((len(L),k)); Cmod=np.zeros((len(L),k,k)); T=Tc[:,0]
        for il in range(len(L)):
            Cmod[il]=Ta[il]+np.diag(shots)
            if not np.all(np.isfinite(Cmod[il])) or Ta[il,0,0]<=0: continue
            weights[il]=np.linalg.solve(Cmod[il],Tc[il])
        kap=combine_maps([m*mask for m in maps],weights,edges,pix)
        combined+=kap
        out['slices'].append({'zmin':s['zmin'],'zmax':s['zmax'],'tracers':names,'map':kap,'mask':mask,'weights':weights,
                              'Chat':Chat,'Cmodel':Cmod,'S_unit_bias':T,'S_theory_L':Lth,'S_theory':S_ll,'S_cross_theory':S_cross,'shot_s':shots,'mock2d':mock2d})
    out['combined']=combined; out['common_mask']=common_mask
    return out


def lowz_bundles(mock,cfg,fixed=None,n=128):
    """Band templates at the sightlines: per slice, combined, truth (kappa_lya), and the fixed (other-realisation)
    slice/combined maps sampled at this mock's sightlines."""
    from templates import flat_sky_band_templates
    sl=mock.sightlines; tm=tracer_maps(mock,cfg,n=n); pix=tm['pixel_size']
    b={'maps':tm,'templates':{}}
    for si,s in enumerate(tm['slices']):
        b['templates'][f'slice{si}'],_=flat_sky_band_templates(s['map'],sl.ra,sl.dec,pix,source=f"slice {s['zmin']}-{s['zmax']} Wiener-combined tracers")
    b['templates']['combined'],_=flat_sky_band_templates(tm['combined'],sl.ra,sl.dec,pix,source='sum of Wiener-filtered slice maps')
    tp=tm['side_rad']/mock.maps['kappa_lya'].shape[0]
    b['templates']['truth'],_=flat_sky_band_templates(mock.maps['kappa_lya'],sl.ra,sl.dec,tp,source='unfiltered truth')
    if fixed is not None:
        b['templates']['fixed_combined'],_=flat_sky_band_templates(fixed['combined'],sl.ra,sl.dec,fixed['pixel_size'],source='combined map of another realisation')
        for si,s in enumerate(fixed['slices']):
            b['templates'][f'fixed_slice{si}'],_=flat_sky_band_templates(s['map'],sl.ra,sl.dec,fixed['pixel_size'],source='slice map of another realisation')
    return b


def save_tracer_maps(path,tm,group='lowz_maps'):
    import h5py
    with h5py.File(path,'a') as f:
        if group in f: del f[group]
        g=f.create_group(group); g.attrs['pixel_size']=tm['pixel_size']; g.attrs['side_rad']=tm['side_rad']
        g['combined']=tm['combined']; g['common_mask']=tm['common_mask']; g['edges']=tm['edges']; g['L']=tm['L']
        for si,s in enumerate(tm['slices']):
            gg=g.create_group(f'slice{si}')
            for k in ('map','mask','weights','Chat','Cmodel','S_unit_bias','S_theory_L','S_theory'): gg[k]=s[k]
            gg.attrs['meta']=json.dumps({k:s[k] for k in ('zmin','zmax','tracers','shot_s')})
        for name,t in tm['tracers'].items():
            gg=g.create_group(f'tracer_{name}')
            for k in ('auto_spectrum','theory_unit_bias'): gg[k]=t[k]
            gg.attrs['meta']=json.dumps({k:v for k,v in t.items() if k not in ('auto_spectrum','theory_unit_bias')},default=float)


def tracer_summary(tm):
    """JSON-able summary of the template construction (bias fits, weights, shot noise)."""
    band=(tm['L']>=SCIENCE_L[0])&(tm['L']<=SCIENCE_L[1])
    def ratio(s):
        d=np.array([np.diag(s['Chat'][il]) for il in np.flatnonzero(band)]); m=np.array([np.diag(s['Cmodel'][il]) for il in np.flatnonzero(band)])
        return (d.sum(axis=0)/np.maximum(m.sum(axis=0),1e-30)).tolist()
    return {'tracers':{name:{k:v for k,v in t.items() if k not in ('auto_spectrum','theory_unit_bias')} for name,t in tm['tracers'].items()},
            'slices':[{'zmin':s['zmin'],'zmax':s['zmax'],'tracers':s['tracers'],'shot_s':s['shot_s'],
                       'weights_science_band':s['weights'][band].tolist(),'measured_over_model_auto_band':ratio(s)} for s in tm['slices']]}


# ----------------------------------------------------------------------------------------------------------------
# Amplitudes: per slice, combined, and the jackknife combination of the slice amplitudes.

def joint_amplitudes(cat,bundle,cfg,sl,names):
    """Fit each named template set on the same regions; return per-name summaries and the joint jackknife
    covariance of the science amplitudes (regions aligned across fits) with the optimal combination."""
    from run_mock_validation import midpoint_regions, common_science
    from amplitude import amplitude
    reg,nside,nreg=midpoint_regions(cat,sl); out={}; jk=[]
    for name in names:
        r=amplitude(cat,bundle['templates'][name],cfg.g1,reg); s=common_science(r)
        out[name]={'A':s['A'],'jk_error':s['jk_error'],'sigma_F':s['sigma_F'],'curl':float(np.mean(r.A[3:6])),
                   'curl_jk_error':float(np.sqrt(np.mean(r.jk_error[3:6]**2))),'nside_jk':nside,'nregion':nreg,'bands':r.A.tolist()}
        jk.append(np.asarray(s['jk']))
    return out,np.asarray(jk),reg


def optimal_combination(A,jk):
    """A: slice amplitudes [k]; jk: jackknife samples [k, nregions]. C = (n-1)/n sum (jk - mean)(jk - mean)^T,
    Hartlap-corrected inverse; A_opt = (1^T C^-1 A)/(1^T C^-1 1)."""
    A=np.asarray(A,float); jk=np.asarray(jk,float); k,nr=jk.shape
    d=jk-jk.mean(axis=1,keepdims=True); C=(nr-1)/nr*d@d.T
    if nr<=k+2: return {'A':float('nan'),'error':float('nan'),'covariance':C.tolist(),'note':'too few regions'}
    hart=(nr-k-2)/(nr-1); Ci=hart*np.linalg.pinv(C); one=np.ones(k)
    var=1/float(one@Ci@one); Aopt=float(one@Ci@A*var)
    naive_w=1/np.maximum(np.diag(C),1e-30); naive=float(np.sum(naive_w*A)/np.sum(naive_w))
    corr=C/np.sqrt(np.outer(np.diag(C),np.diag(C)))
    return {'A':Aopt,'error':float(np.sqrt(var)),'naive_A':naive,'naive_error':float(1/np.sqrt(np.sum(naive_w))),
            'covariance':C.tolist(),'correlation':corr.tolist(),'hartlap':hart,'nregions':nr}
