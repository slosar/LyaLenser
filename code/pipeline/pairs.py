"""Sightline-pair discovery and deterministic 11-by-6 compression kernel."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import argparse, json, time
import numpy as np
from scipy.spatial import cKDTree
from numba import njit, prange, get_num_threads, set_num_threads

try:
    from .config import Config, SightlineSet
except ImportError:
    from config import Config, SightlineSet

ACCUMULATORS = ("v0", "v1", "vc", "m0", "m1", "m2", "mc", "mcc",
                "beta0", "beta1", "betac")


@dataclass
class PairCatalogue:
    a: np.ndarray
    b: np.ndarray
    thx: np.ndarray
    thy: np.ndarray
    theta: np.ndarray
    accum: np.ndarray
    npair: np.ndarray
    attrs: dict = field(default_factory=dict)

    def __post_init__(self):
        self.a = np.asarray(self.a, np.int32); self.b = np.asarray(self.b, np.int32)
        self.thx = np.asarray(self.thx, np.float32); self.thy = np.asarray(self.thy, np.float32)
        self.theta = np.asarray(self.theta, np.float32)
        self.accum = np.asarray(self.accum)
        self.npair = np.asarray(self.npair, np.int32)
        if self.accum.ndim != 3 or self.accum.shape[1:] != (11, 6):
            raise ValueError("accum must be [Npair,11,6]")

    def named(self, name):
        return self.accum[:, ACCUMULATORS.index(name), :]

    def save(self, path, group="all"):
        import h5py
        path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "a") as f:
            if group in f: del f[group]
            g = f.create_group(group)
            for k in ("a","b","thx","thy","theta","npair"):
                g.create_dataset(k, data=getattr(self, k), compression="gzip", shuffle=True)
            g.create_dataset("accum", data=self.accum.astype(np.float32), compression="gzip", shuffle=True)
            g.attrs["accumulator_order"] = json.dumps(ACCUMULATORS)
            g.attrs["accumulation_precision"] = "float64; storage float32"
            for k,v in self.attrs.items():
                if np.isscalar(v) or isinstance(v, str): g.attrs[k] = v


def _unit_vectors(ra, dec):
    r = np.deg2rad(ra); d = np.deg2rad(dec); c = np.cos(d)
    return np.column_stack((c*np.cos(r), c*np.sin(r), np.sin(d)))


def pair_geometry(ra, dec, a, b):
    """Local direction theta_ab = theta_a-theta_b, from b to a."""
    dra = (np.deg2rad(ra[a]-ra[b]) + np.pi) % (2*np.pi) - np.pi
    dmid = .5*np.deg2rad(dec[a]+dec[b])
    dx = dra*np.cos(dmid)
    dy = np.deg2rad(dec[a]-dec[b])
    theta = np.hypot(dx,dy)
    safe = np.where(theta > 0, theta, 1)
    return (dx/safe).astype(np.float32), (dy/safe).astype(np.float32), theta.astype(np.float32)


def find_pairs(sl: SightlineSet, theta_max: float):
    tree = cKDTree(_unit_vectors(sl.ra, sl.dec))
    ij = tree.query_pairs(2*np.sin(theta_max/2), output_type="ndarray")
    if len(ij) == 0:
        e = np.empty(0, np.int32); f = np.empty(0, np.float32)
        return e,e.copy(),f,f.copy(),f.copy()
    a = ij[:,0].astype(np.int32); b = ij[:,1].astype(np.int32)
    thx,thy,theta = pair_geometry(sl.ra,sl.dec,a,b)
    keep = theta <= theta_max
    return a[keep],b[keep],thx[keep],thy[keep],theta[keep]


@njit(cache=True, inline="always")
def _interp(rp, rz, rp0, drp, nrp, rz0, drz, nrz, xi, xirp):
    if rp < rp0 or rz < rz0 or rp > rp0+drp*(nrp-1) or rz > rz0+drz*(nrz-1):
        return 0.0,0.0
    i=min(int((rp-rp0)/drp),nrp-2); j=min(int((rz-rz0)/drz),nrz-2)
    x=(rp-(rp0+i*drp))/drp; y=(rz-(rz0+j*drz))/drz
    k=i*nrz+j; k2=(i+1)*nrz+j
    v=((1-x)*(1-y)*xi[k]+x*(1-y)*xi[k2]+(1-x)*y*xi[k+1]+x*y*xi[k2+1])
    g=((1-x)*(1-y)*xirp[k]+x*(1-y)*xirp[k2]+(1-x)*y*xirp[k+1]+x*y*xirp[k2+1])
    return v,g


@njit(cache=True, inline="always")
def _interp_layer(rp, rz, cm, rp0, drp, nrp, rz0, drz, nrz, chi0, dchi, nchi, xi, xirp):
    """`_interp` on a layered table (xi, xirp flattened [nchi, nrp, nrz]): linear between the two chi layers that
    bracket the pair's mean distance ``cm``, clamped outside the node range; nchi == 1 is the plain 2-D table."""
    n2 = nrp*nrz
    if nchi <= 1:
        return _interp(rp, rz, rp0, drp, nrp, rz0, drz, nrz, xi, xirp)
    t = (cm-chi0)/dchi
    if t <= 0.0:
        return _interp(rp, rz, rp0, drp, nrp, rz0, drz, nrz, xi[:n2], xirp[:n2])
    if t >= nchi-1:
        return _interp(rp, rz, rp0, drp, nrp, rz0, drz, nrz, xi[(nchi-1)*n2:], xirp[(nchi-1)*n2:])
    i = int(t); f = t-i
    v0, g0 = _interp(rp, rz, rp0, drp, nrp, rz0, drz, nrz, xi[i*n2:(i+1)*n2], xirp[i*n2:(i+1)*n2])
    v1, g1 = _interp(rp, rz, rp0, drp, nrp, rz0, drz, nrz, xi[(i+1)*n2:(i+2)*n2], xirp[(i+1)*n2:(i+2)*n2])
    return (1.0-f)*v0+f*v1, (1.0-f)*g0+f*g1


@njit(cache=True, inline="always")
def _shape_bin(rp, rz):
    if rp < 0 or rp > 30 or rz < 0 or rz > 30: return -1
    ir = 0 if rp < 10 else (1 if rp < 20 else 2)
    iz = 0 if rz < 10 else 1
    return 2*ir+iz


@njit(parallel=True, cache=True)
def _accumulate_kernel(pix_start, chi, delta, weight, pa, pb, theta,
                       rp_grid, rz_grid, xi, xirp, rpmax, rzmax, chi_ref, slab, slab_edges, slab_index, rpmin,
                       theta_true, expectation, region, chi0, dchi, nchi):
    """Pair sums. With ``expectation`` the product delta_p delta_q is replaced by the table's xi at the TRUE
    separation (``theta_true`` per sightline pair, the unshifted geometry of an injection) while the kernel, the
    selection and the mean field use the observed (shifted) geometry: the noise-free expectation of an injection."""
    n=pa.size; out=np.zeros((n,11,6),np.float64); counts=np.zeros(n,np.int32)
    rp0=rp_grid[0]; rz0=rz_grid[0]
    drp=rp_grid[1]-rp_grid[0]; drz=rz_grid[1]-rz_grid[0]
    nrp=rp_grid.size; nrz=rz_grid.size
    for ip in prange(n):
        a=pa[ip]; b=pb[ip]
        qb=pix_start[b]; qend=pix_start[b+1]
        for p in range(pix_start[a],pix_start[a+1]):
            cp=np.float64(chi[p])
            while qb < qend and np.float64(chi[qb]) < cp-rzmax: qb += 1
            q=qb
            while q < qend and np.float64(chi[q]) <= cp+rzmax:
                cq=np.float64(chi[q]); dc=cp-cq; rz=abs(dc); cm=.5*(cp+cq)
                rp=cm*np.float64(theta[ip])
                selected = slab[p]>=0 and slab[q]>=0
                if region[p]>0 and region[q]>0: selected = False   # B x B carries Lyb-Lyb, not modelled
                if slab_index>=0:
                    selected = selected and slab_edges[slab_index,0]<=cm<slab_edges[slab_index,1]
                if rp <= rpmax and rp >= rpmin and selected:
                    ib=_shape_bin(rp,rz)
                    if ib >= 0:
                        xv,xg=_interp_layer(rp,rz,cm,rp0,drp,nrp,rz0,drz,nrz,chi0,dchi,nchi,xi,xirp)
                        G=cm*xg; dm=cm-chi_ref
                        ww=np.float64(weight[p])*np.float64(weight[q])
                        if expectation:
                            dd,_=_interp_layer(cm*np.float64(theta_true[ip]),rz,cm,rp0,drp,nrp,rz0,drz,nrz,chi0,dchi,nchi,xi,xirp)
                        else:
                            dd=np.float64(delta[p])*np.float64(delta[q])
                        out[ip,0,ib]+=ww*dd*G
                        out[ip,1,ib]+=ww*dd*G*dm
                        out[ip,2,ib]+=ww*dd*G*dc*.5
                        gg=ww*G*G
                        out[ip,3,ib]+=gg; out[ip,4,ib]+=gg*dm; out[ip,5,ib]+=gg*dm*dm
                        out[ip,6,ib]+=gg*dc; out[ip,7,ib]+=gg*dc*dc*.25
                        bx=ww*xv*G
                        out[ip,8,ib]+=bx; out[ip,9,ib]+=bx*dm; out[ip,10,ib]+=bx*dc*.5
                        counts[ip]+=1
                q+=1
    return out,counts


def accumulate(sl, pairs, xi_table, cfg: Config, shifted_positions=None, true_positions=None):
    """Pair sums on the observed geometry (``shifted_positions`` (ra, dec) per sightline replaces sl.ra/dec).
    ``true_positions`` switches to the expectation mode: delta_p delta_q -> xi(true separation) while the kernel,
    the selection and the mean field use the observed geometry (noise-free injection expectation)."""
    a,b,thx,thy,theta = pairs
    if shifted_positions is not None:
        pos=np.asarray(shifted_positions)
        thx,thy,theta=pair_geometry(pos[:,0],pos[:,1],a,b)
    expectation=true_positions is not None
    if expectation:
        pos=np.asarray(true_positions); _,_,theta_true=pair_geometry(pos[:,0],pos[:,1],a,b)
    else: theta_true=theta
    region=np.asarray(getattr(sl,"region",None) if getattr(sl,"region",None) is not None
                      else np.zeros(len(sl.chi),np.int8),np.int8)
    chi0,dchi,nchi=xi_table.layers() if hasattr(xi_table,"layers") else (0.,1.,1)
    out,n=_accumulate_kernel(sl.pix_start,sl.chi,sl.delta,sl.w,a,b,theta,
                             xi_table.r_perp,xi_table.r_par,xi_table.xi.ravel(),
                             xi_table.xi_rp.ravel(),cfg.r_perp_max,cfg.r_par_max,cfg.chi_ref,sl.slab,
                             np.asarray([[float(__import__("cosmo").chi(z)) for z in bounds] for bounds in cfg.slabs]),cfg.slab_index,
                             float(getattr(cfg,"r_perp_min",0.0)),theta_true,expectation,region,chi0,dchi,nchi)
    keep=n>0
    attrs={"chi_ref":cfg.chi_ref,"accumulation_precision":"float64","region_B_pixels":int((region>0).sum()),
           "xi_layers":int(nchi),
           "storage_precision":"float32","pair_direction":"theta_a-theta_b","r_perp_min":float(getattr(cfg,"r_perp_min",0.0)),
           "expectation":expectation}
    return PairCatalogue(a[keep],b[keep],thx[keep],thy[keep],theta[keep],out[keep],n[keep],attrs)


def pair_midpoint_regions(cat, sl, nside=8):
    import healpy as hp
    va=_unit_vectors(sl.ra[cat.a],sl.dec[cat.a]); vb=_unit_vectors(sl.ra[cat.b],sl.dec[cat.b])
    v=va+vb; v/=np.linalg.norm(v,axis=1)[:,None]
    return hp.vec2pix(nside,v[:,0],v[:,1],v[:,2]).astype(np.int32)


def benchmark(sl, fraction=1.0, xi_table=None, cfg=None, output=None):
    cfg=Config() if cfg is None else cfg
    if xi_table is None:
        from xi_model import XiTable
        g=np.arange(0,40.25,.25); x=np.exp(-np.hypot(g[:,None],g[None,:])/10)
        xi_table=XiTable(g,g,x,np.gradient(x,.25,axis=0))
    theta_max=cfg.r_perp_max/max(float(sl.chi.min()),1)
    pairs=find_pairs(sl,theta_max)
    m=max(1,int(len(pairs[0])*fraction)) if len(pairs[0]) else 0
    pairs=tuple(x[:m] for x in pairs)
    accumulate(sl,pairs,xi_table,cfg) # compile/warm
    t=time.perf_counter(); cat=accumulate(sl,pairs,xi_table,cfg); wall=time.perf_counter()-t
    pixels=int(cat.npair.sum()); threads=get_num_threads()
    result={"sightline_pairs":len(cat.a),"pixel_pairs":pixels,"threads":threads,"wall_s":wall,
            "pixel_pairs_per_s":pixels/max(wall,1e-12),
            "pixel_pairs_per_s_per_thread":pixels/max(wall*threads,1e-12),
            "extrapolated_full_wall_s":wall/max(fraction,1e-12)}
    if output:
        Path(output).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2)); return result


if __name__ == "__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("mock"); ap.add_argument("--fraction",type=float,default=.1)
    args=ap.parse_args()
    from mock import load_sightlines
    benchmark(load_sightlines(args.mock),args.fraction)


def accumulate_slabs(sl,pairs,xi_table,cfg,path=None):
    cats={}
    for i in range(len(cfg.slabs)):
        cat=accumulate(sl,pairs,xi_table,cfg.copy(slab_index=i))
        cats[f"slab{i}"]=cat
        if path is not None: cat.save(path,f"slab{i}")
    return cats
