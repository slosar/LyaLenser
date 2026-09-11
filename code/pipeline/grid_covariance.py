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


def xi_from_mock_grid(cfg,shape,dx=2.,dz=.5,nangle=128,return_cube=False):
    """Same discrete P(k) and Fourier normalization as mock._fft_fields."""
    import gc
    nx,ny,nz=map(int,shape); pf=ForestPower(model='kaiser'); p=pf.p
    kt=np.hypot(2*np.pi*fftfreq(nx,dx)[:,None],2*np.pi*fftfreq(ny,dx)[None,:])
    kz=2*np.pi*rfftfreq(nz,dz); kg=np.geomspace(1e-4,49,3000); pg=pf.plin(kg)
    spec=np.empty((nx,ny,len(kz)),np.complex64)
    for iz,k in enumerate(kz):
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
