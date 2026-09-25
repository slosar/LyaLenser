"""Expected covariance of the mock grid after periodic trilinear sampling.

C_grid[n] = IFFT(P_F(k) / cell_volume)[n]. For uniform subcell phase,
C_interp(r) = sum_n C_grid[n] product_i B3(r_i/d_i - n_i), where B3 is
Lambda * Lambda, the autocorrelation of the linear interpolation hat.
This is the exact discrete-grid/pixel-window calculation, including aliases;
no empirical scalar is applied. A uniform transverse-angle average supplies
the same two-coordinate covariance interface as the production xi table.
"""
import numpy as np
from numba import njit,prange
from scipy.fft import fftfreq,rfftfreq,irfftn
from forest_power import ForestPower
from xi_model import XiTable

@njit(cache=True)
def b3(x):
    a=abs(x); sign=1. if x>=0 else -1.
    if a<1: return (4-6*a*a+3*a*a*a)/6,sign*(-2*a+1.5*a*a)
    if a<2: return (2-a)**3/6,-.5*sign*(2-a)**2
    return 0.,0.

@njit(parallel=True,cache=True)
def sample_covariance(cube,grid,dx,dz,nangle):
    n=len(grid); xi=np.zeros((n,n)); grad=np.zeros((n,n))
    hx,hy,hz=cube.shape[0]//2,cube.shape[1]//2,cube.shape[2]//2
    for ir in prange(n):
        for iz in range(n):
            z=grid[iz]/dz; kz=int(np.floor(z)); value=0.; derivative=0.
            for angle in range(nangle):
                theta=2*np.pi*(angle+.5)/nangle; ct=np.cos(theta); st=np.sin(theta)
                x=grid[ir]*ct/dx; y=grid[ir]*st/dx; ix=int(np.floor(x)); iy=int(np.floor(y))
                for i in range(ix-1,ix+3):
                    bx,dbx=b3(x-i)
                    for j in range(iy-1,iy+3):
                        by,dby=b3(y-j)
                        for k in range(kz-1,kz+3):
                            bz,_=b3(z-k); cov=cube[i+hx,j+hy,k+hz]
                            value+=cov*bx*by*bz
                            derivative+=cov*(dbx*ct*by+bx*dby*st)*bz/dx
            xi[ir,iz]=value/nangle; grad[ir,iz]=derivative/nangle
    return xi,grad


def xi_from_mock_grid(cfg,shape,dx=2.,dz=.5,nangle=128,return_cube=False,power=None):
    """Same discrete P(k) and Fourier normalization as mock._fft_fields.

    ``power(k_par, k_perp)`` replaces the fiducial Kaiser forest spectrum (used for the basis spectra of
    ``xi_fit``); it must accept broadcasting arrays.
    """
    import gc
    nx,ny,nz=map(int,shape); pf=ForestPower(model='kaiser'); p=pf.p
    kt=np.hypot(2*np.pi*fftfreq(nx,dx)[:,None],2*np.pi*fftfreq(ny,dx)[None,:])
    kz=2*np.pi*rfftfreq(nz,dz); kg=np.geomspace(1e-4,49,3000); pg=pf.plin(kg)
    spec=np.empty((nx,ny,len(kz)),np.complex64)
    for iz,k in enumerate(kz):
        if power is not None:
            spec[:,:,iz]=power(k,kt)/(dx*dx*dz); continue
        kk=np.hypot(kt,k); pl=np.interp(np.clip(kk,kg[0],kg[-1]),kg,pg)
        mu2=k*k/np.maximum(kk*kk,1e-30)
        response=p['b_F']*(1+p['beta_F']*mu2)*np.exp(-.5*(kk/p['kp'])**2)
        spec[:,:,iz]=pl*response**2/(dx*dx*dz)
    spec[0,0,0]=0
    covariance=irfftn(spec,s=(nx,ny,nz),workers=4,overwrite_x=True); del spec
    hx=int(np.ceil(cfg.xi_max/dx))+3; hz=int(np.ceil(cfg.xi_max/dz))+3
    ix=np.arange(-hx,hx+1)%nx; iy=np.arange(-hx,hx+1)%ny; iz=np.arange(-hz,hz+1)%nz
    cube=np.ascontiguousarray(covariance[np.ix_(ix,iy,iz)],dtype=np.float64)
    del covariance; gc.collect()
    grid=np.arange(0,cfg.xi_max+cfg.xi_step*.5,cfg.xi_step)
    xi,grad=sample_covariance(cube,grid,dx,dz,nangle)
    table=XiTable(grid,grid,xi,grad,{'provider':'discrete_mock_grid','shape':list(shape),'dx':dx,'dz':dz,'angular_samples':nangle,
                                    'pixel_window':'exact phase-averaged product of cubic B3 autocorrelations; aliases retained',
                                    'empirical_factor':1.0})

    return (table,cube) if return_cube else table


def projected_grid_table(raw, sl, cfg):
    """Expected measured table for uniform-weight mock forests, including P C P^T.

    The mock has the same radial grid on every forest and constant weight
    within a forest. This evaluates its continuum projection before binning
    and before the identical measured-table smoothing. The transverse angle
    average remains the approximation documented for the raw grid table.
    """
    from scipy.interpolate import RegularGridInterpolator
    from xi_model import xi_from_counts
    c=np.asarray(sl.chi[sl.pix_start[0]:sl.pix_start[1]],float)
    rz=abs(c[:,None]-c[None,:]); u=np.column_stack((np.ones(len(c)),c-c.mean()))
    v=u@np.linalg.inv(u.T@u)
    bins=rz.astype(int); valid=bins<cfg.xi_max
    counts=np.bincount(bins[valid],minlength=int(cfg.xi_max))
    fn=RegularGridInterpolator((raw.r_perp,raw.r_par),raw.xi,bounds_error=False,fill_value=0.)
    num=np.zeros((int(cfg.xi_max),int(cfg.xi_max))); den=np.zeros_like(num)
    # Average across each transverse bin as well as the actual radial pixels.
    for ir in range(len(num)):
        row=np.zeros(len(counts))
        for rp in (ir+.125,ir+.375,ir+.625,ir+.875):
            matrix=fn((np.full_like(rz,rp),rz))
            mv=matrix@v
            projected=matrix-u@(v.T@matrix)-mv@u.T+u@(v.T@mv)@u.T
            row+=np.bincount(bins[valid],weights=projected[valid],minlength=len(counts))/4
        num[ir]=row; den[ir]=counts
    return xi_from_counts(num,den,cfg)


def derivative_comparison(measured,predicted):
    """Response-integral projection; geometric r_perp^3 weight, r_parallel <=30.

    Report a signed normalization residual and an orthogonal shape residual;
    neither uses paired lensing amplitudes or fits a correction to validation.
    """
    rp=predicted.r_perp; rz=predicted.r_par
    use=(rp>5)&(rp<30); radial=rz<=30
    w=rp[use,None]**3
    p=predicted.xi_rp[use][:,radial]; m=measured.xi_rp[use][:,radial]
    norm=np.sum(w*p*p); coefficient=float(np.sum(w*m*p)/norm)
    return {'coefficient':coefficient,'residual':coefficient-1,
            'relative_norm':float(np.sqrt(np.sum(w*(m-p)**2)/norm)),
            'range':'5 < r_perp < 30; 0 <= r_parallel <= 30',
            'weight':'r_perp^3 times predicted derivative squared'}
