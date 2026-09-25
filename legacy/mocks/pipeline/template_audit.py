"""Map-level matched-template audit with equal radial and angular operators.

Native maps are samples at grid nodes. Catalogue objects represent cells
CENTERED on those nodes. Exact overlap integration, including periodic edges,
places these cells on the count map; Fourier resampling is a different operator.
"""
import numpy as np
from scipy.ndimage import zoom

MARGINS = (0, 150, 300)

def overlap_matrix(nin, nout, offset=0.):
    # Native grid nodes are at i/nin; output pixel edges are at j/nout.
    left=(np.arange(nin)-.5+offset)/nin
    right=left+1/nin
    lo=np.arange(nout)/nout; hi=lo+1/nout
    w=np.zeros((nout,nin))
    for shift in (-1,0,1):
        w+=np.maximum(0,np.minimum(hi[:,None],right[None,:]+shift)-
                      np.maximum(lo[:,None],left[None,:]+shift))*nout
    return w

def count_pixel_map(a,n=128,offset=0.):
    wx=overlap_matrix(a.shape[0],n,offset); wy=overlap_matrix(a.shape[1],n,offset)
    return wx@np.asarray(a,float)@wy.T

def continuous_projections(density,chis,kernel,cforest,correction,bias,tracer=None):
    """Persist deterministic intensity and density integrals, without Poisson noise.

    The finite-volume lognormal realization is normalized radially just as the
    catalogue nbar estimator is. Both nominal and measured-mean versions are
    retained to isolate that integral constraint from cell centering.
    ``tracer`` is the (radially smoothed) field the lognormal is applied to;
    the truth integrals always use ``density``.
    """
    tracer=density if tracer is None else tracer
    dz=float(chis[1]-chis[0]); maps={}
    for margin in MARGINS:
        lower,upper=cforest[0]-margin,cforest[1]+margin
        weight=np.maximum(0,np.minimum(chis+dz/2,upper)-np.maximum(chis-dz/2,lower))*kernel
        truth=np.zeros(density.shape[:2]); intensity=np.zeros_like(truth); normalized=np.zeros_like(truth)
        indices=np.flatnonzero(weight)
        edges=np.linspace(lower,upper,41)
        radial_bin=np.clip(np.searchsorted(edges,chis,side='right')-1,0,39)
        sums=np.zeros(40); volumes=np.zeros(40)
        for i in indices:
            lam=np.exp(np.clip(bias*np.asarray(tracer[:,:,i],float)-correction,-30,30))
            width=weight[i]/kernel[i]
            sums[radial_bin[i]]+=lam.mean()*width
            volumes[radial_bin[i]]+=width
        means=np.divide(sums,volumes,out=np.ones_like(sums),where=volumes>0)
        for i in indices:
            plane=np.asarray(density[:,:,i],float)
            lam=np.exp(np.clip(bias*np.asarray(tracer[:,:,i],float)-correction,-30,30))
            truth+=weight[i]*plane
            intensity+=weight[i]*(lam-1)/bias
            normalized+=weight[i]*(lam/means[radial_bin[i]]-1)/bias
        maps[f'kappa_range_{margin}']=truth.astype('f4')
        maps[f'intensity_range_{margin}']=intensity.astype('f4')
        maps[f'intensity_normalized_range_{margin}']=normalized.astype('f4')
    return maps

def regression(template,truth,mask):
    use=np.asarray(mask)>0
    x=np.asarray(truth)[use]; y=np.asarray(template)[use]
    # Free intercept removes the unobservable angular monopole.
    x=x-x.mean(); y=y-y.mean()
    if len(x)<3 or np.dot(x,x)<=0: raise ValueError('empty or constant template audit truth')
    return float(np.dot(x,y)/np.dot(x,x))


def band_regression(template,truth,mask,side_rad,lmin=40.,lmax=300.):
    """Cross/auto coefficient of the masked maps restricted to the science band.

    The full-resolution pixel regression is dominated by the smallest scales;
    the deprojection only uses 40 <= L <= 300, so the gate is judged there.
    """
    m=np.asarray(mask,float); n=m.shape[0]
    t=np.fft.fft2((np.asarray(template,float)-np.mean(template[m>0]))*m)
    k=np.fft.fft2((np.asarray(truth,float)-np.mean(truth[m>0]))*m)
    L=2*np.pi*np.fft.fftfreq(n,d=side_rad/n); LL=np.hypot(L[:,None],L[None,:])
    sel=(LL>=lmin)&(LL<=lmax)
    if sel.sum()<2 or np.sum(abs(k[sel])**2)<=0: raise ValueError('empty science band in template audit')
    return float(np.sum((t*np.conj(k)).real[sel])/np.sum(abs(k[sel])**2))

def audit_mock(mock,path,bundle_factory):
    import h5py
    result={}
    with h5py.File(path,'a') as f:
        if 'template_audit' in f: del f['template_audit']
        group=f.create_group('template_audit')
        for margin in MARGINS:
            b=bundle_factory(mock,template_margin=margin)
            mask=b['common_mask']; n=mask.shape[0]
            truth=count_pixel_map(mock.maps[f'kappa_range_{margin}'],n)
            nominal=count_pixel_map(mock.maps[f'intensity_range_{margin}'],n)
            realspace=count_pixel_map(mock.maps[f'intensity_normalized_range_{margin}'],n)
            continuous=count_pixel_map(mock.maps[f'continuous_catalogue_range_{margin}'],n)
            old=count_pixel_map(mock.maps[f'intensity_normalized_range_{margin}'],n,offset=.5)
            # Matched output already contains the common mask; compare identical operators.
            groupm=group.create_group(str(margin))
            for k,v in dict(truth=truth,continuous=continuous,realspace=realspace,nominal=nominal,old_offset=old,
                            sampled=b['matched_map'],mask=mask).items(): groupm[k]=v
            row={k:regression(v*mask,truth*mask,mask) for k,v in
                 [('continuous',continuous),('realspace',realspace),('nominal',nominal),('old_offset',old)]}
            row['sampled']=regression(b['matched_map'],truth*mask,mask)
            # band_regression applies the common mask itself, exactly once, to an UNMASKED map: the sampled
            # template therefore enters unmasked (review 5, finding 5: the masked map was masked twice).
            side=float(b['side_rad'])
            for k,v in [('continuous',continuous),('realspace',realspace),('nominal',nominal),('old_offset',old),('sampled',b['matched_map_unmasked'])]:
                row[k+'_band']=band_regression(v,truth,mask,side)
            row['chi_range']=[float(__import__('cosmo').chi(mock.sightlines.attrs['zmin']))-margin,
                              float(__import__('cosmo').chi(mock.sightlines.attrs['zmax']))+margin]
            result[str(margin)]=row
    return result


def continuous_catalogue_projection(density,rsd,chis,dx,dz,cref,cforest,correction,
                                    completeness,magnification_map,bias=3.5):
    """No-Poisson catalogue expectation, including radial RSD and actual selection.

    ``density`` is the field the lognormal intensity is built from (the radially
    smoothed tracer field in the campaign mocks). Radial transport is evaluated
    at native cell centers; boundary cells use fractional overlap. Catalogue
    nbar uses the same 40 observed-chi bins. Exact angular cell overlap is
    applied later, identically to range truth.
    """
    from mock import ray_points,_interp3_chunked
    from cross_spectrum import kernel,Z_CMB
    from cosmo import chi as chi_of_z
    nx,ny,nz=density.shape; size=nx*ny
    ra=180+np.rad2deg((np.arange(nx)-nx/2)*dx/cref)/np.cos(np.deg2rad(30))
    dec=30+np.rad2deg((np.arange(ny)-ny/2)*dx/cref)
    rr,dd=np.meshgrid(ra,dec,indexing='ij')
    comp=np.asarray(completeness,float).ravel(); compmean=comp.mean()
    mag=np.ones(size) if magnification_map is None else np.clip(1+.5*np.asarray(magnification_map).ravel(),.05,None)
    totals={m:np.zeros(40) for m in MARGINS}
    numerators={m:np.zeros((40,size)) for m in MARGINS}
    cmbchi=float(chi_of_z(Z_CMB)); angular=np.arange(size)
    radial_kernel=kernel(chis,cmbchi)
    for start in range(0,nz,16):
        stop=min(start+16,nz); cs=chis[start:stop]
        points=ray_points(rr.ravel()[:,None],dd.ravel()[:,None],cs[None,:],dx,dz,rsd.shape,chis[0])
        shift=_interp3_chunked(rsd,points.reshape(-1,3)).reshape(size,len(cs))
        observed=cs[None,:]+shift
        lam=np.exp(np.clip(bias*np.asarray(density[:,:,start:stop],float).reshape(size,len(cs))-correction,-30,30))*mag[:,None]
        weights=np.interp(observed.ravel(),chis,radial_kernel).reshape(observed.shape)
        for margin in MARGINS:
            lo,hi=cforest[0]-margin,cforest[1]+margin; step=(hi-lo)/40
            bins=np.clip(((observed-lo)/step).astype(int),0,39)
            length=np.maximum(0,np.minimum(observed+dz/2,hi)-np.maximum(observed-dz/2,lo))
            mass=lam*length
            for j in range(len(cs)):
                totals[margin]+=np.bincount(bins[:,j],weights=mass[:,j]*comp,minlength=40)/size
                np.add.at(numerators[margin],(bins[:,j],angular),mass[:,j]*weights[:,j])
    maps={}
    for margin in MARGINS:
        lo,hi=cforest[0]-margin,cforest[1]+margin; step=(hi-lo)/40
        centers=lo+(np.arange(40)+.5)*step
        h=totals[margin]
        if np.any(h<=0): raise ValueError('empty continuous radial selection bin')
        q=compmean/bias*np.sum(numerators[margin]*(step/h)[:,None],axis=0)
        # Same global Nq/Nr random subtraction, with a smooth angular selection.
        random=h.sum()/(bias*(hi-lo))*np.sum(step**2*kernel(centers,cmbchi)/h)
        maps[f'continuous_catalogue_range_{margin}']=(q-random).reshape(nx,ny).astype('f4')
    return maps
