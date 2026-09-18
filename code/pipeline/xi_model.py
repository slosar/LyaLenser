"""Observed forest correlation tables used by the pair response."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
import numpy as np
from scipy.special import j0, j1
from scipy.ndimage import gaussian_filter
from numba import njit, prange

HERE = Path(__file__).resolve().parent
CODE = HERE.parent
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
from forest_power import ForestPower
try:
    from .config import Config, SightlineSet
except ImportError:
    from config import Config, SightlineSet


@dataclass
class XiTable:
    """Correlation table xi(r_perp, r_par) and its r_perp derivative.

    Either a single 2-D table (``chi_nodes`` None, arrays [n_rp, n_rz]) or, with redshift evolution (iteration 10),
    a LAYERED table: ``chi_nodes`` is a uniform grid of comoving distances and the arrays are [n_chi, n_rp, n_rz];
    the pair kernels interpolate linearly between the two layers that bracket the pair's mean distance and clamp
    outside the node range.
    """
    r_perp: np.ndarray
    r_par: np.ndarray
    xi: np.ndarray
    xi_rp: np.ndarray
    meta: dict | None = None
    counts: tuple | None = None
    chi_nodes: np.ndarray | None = None

    def __post_init__(self):
        self.r_perp = np.ascontiguousarray(self.r_perp, dtype=np.float64)
        self.r_par = np.ascontiguousarray(self.r_par, dtype=np.float64)
        self.xi = np.ascontiguousarray(self.xi, dtype=np.float64)
        self.xi_rp = np.ascontiguousarray(self.xi_rp, dtype=np.float64)
        expected = (len(self.r_perp), len(self.r_par))
        if self.chi_nodes is not None:
            self.chi_nodes = np.ascontiguousarray(self.chi_nodes, dtype=np.float64)
            if len(self.chi_nodes) > 1 and not np.allclose(np.diff(self.chi_nodes), self.chi_nodes[1]-self.chi_nodes[0]):
                raise ValueError("chi_nodes must be uniformly spaced")
            expected = (len(self.chi_nodes),) + expected
        if self.xi.shape != expected or self.xi_rp.shape != expected:
            raise ValueError(f"xi arrays must have shape {expected}")

    # ---- layer description for the numba kernels: (chi0, dchi, nchi); a 2-D table is one layer everywhere
    @property
    def layered(self):
        return self.chi_nodes is not None

    def layers(self):
        if self.chi_nodes is None: return 0.0, 1.0, 1
        n = len(self.chi_nodes)
        return float(self.chi_nodes[0]), (float(self.chi_nodes[1]-self.chi_nodes[0]) if n > 1 else 1.0), n

    def at_chi(self, chi):
        """2-D table at one distance (linear interpolation between layers, clamped)."""
        if self.chi_nodes is None: return self
        c0, dc, n = self.layers()
        t = np.clip((float(chi)-c0)/dc, 0, n-1); i = min(int(t), max(n-2, 0)); f = t-i if n > 1 else 0.
        j = min(i+1, n-1)
        return XiTable(self.r_perp, self.r_par, (1-f)*self.xi[i]+f*self.xi[j], (1-f)*self.xi_rp[i]+f*self.xi_rp[j],
                       dict(self.meta or {}, layer_chi=float(chi)))

    def save(self,path,group="xi"):
        import h5py,json
        with h5py.File(path,"a") as f:
            if group in f: del f[group]
            g=f.create_group(group)
            for k in ("r_perp","r_par","xi","xi_rp"): g[k]=getattr(self,k)
            if self.chi_nodes is not None: g["chi_nodes"]=self.chi_nodes
            g.attrs["meta"]=json.dumps(self.meta)
            if self.counts is not None:
                g["coarse_num"]=self.counts[0]; g["coarse_den"]=self.counts[1]

    def interp(self, rp, rz, chi=None):
        if self.chi_nodes is None:
            return bilinear(rp, rz, self.r_perp, self.r_par, self.xi.ravel(), self.xi_rp.ravel())
        if chi is None: raise ValueError("layered table: interp needs the pair mean distance chi")
        c0, dc, n = self.layers()
        return layered_bilinear(rp, rz, float(chi), self.r_perp, self.r_par, c0, dc, n, self.xi.ravel(), self.xi_rp.ravel())


@njit(cache=True)
def bilinear(rp, rz, rp_grid, rz_grid, xi_flat, xirp_flat):
    """Numba-callable interpolation returning (xi, dxi/dr_perp)."""
    nrz = rz_grid.size
    if rp < rp_grid[0] or rz < rz_grid[0] or rp > rp_grid[-1] or rz > rz_grid[-1]:
        return 0.0, 0.0
    drp = rp_grid[1] - rp_grid[0]
    drz = rz_grid[1] - rz_grid[0]
    i = min(int((rp-rp_grid[0])/drp), rp_grid.size-2)
    j = min(int((rz-rz_grid[0])/drz), rz_grid.size-2)
    tx = (rp-rp_grid[i])/drp
    ty = (rz-rz_grid[j])/drz
    k00 = i*nrz+j; k10 = (i+1)*nrz+j
    def one(a):
        return ((1-tx)*(1-ty)*a[k00] + tx*(1-ty)*a[k10] +
                (1-tx)*ty*a[k00+1] + tx*ty*a[k10+1])
    return one(xi_flat), one(xirp_flat)


@njit(cache=True)
def layered_bilinear(rp, rz, chi, rp_grid, rz_grid, chi0, dchi, nchi, xi_flat, xirp_flat):
    """(xi, dxi/dr_perp) of a layered table at the pair mean distance ``chi``: bilinear in (r_perp, r_par) on the
    two bracketing chi layers, linear between them, clamped to the node range. nchi == 1 is the plain 2-D table."""
    if nchi <= 1:
        return bilinear(rp, rz, rp_grid, rz_grid, xi_flat, xirp_flat)
    n2 = rp_grid.size*rz_grid.size
    t = (chi-chi0)/dchi
    if t <= 0.0:
        return bilinear(rp, rz, rp_grid, rz_grid, xi_flat[:n2], xirp_flat[:n2])
    if t >= nchi-1:
        return bilinear(rp, rz, rp_grid, rz_grid, xi_flat[(nchi-1)*n2:], xirp_flat[(nchi-1)*n2:])
    i = int(t); f = t-i
    v0, g0 = bilinear(rp, rz, rp_grid, rz_grid, xi_flat[i*n2:(i+1)*n2], xirp_flat[i*n2:(i+1)*n2])
    v1, g1 = bilinear(rp, rz, rp_grid, rz_grid, xi_flat[(i+1)*n2:(i+2)*n2], xirp_flat[(i+1)*n2:(i+2)*n2])
    return (1-f)*v0+f*v1, (1-f)*g0+f*g1


def _trap_weights_log(k):
    x = np.log(k)
    w = np.empty_like(x)
    w[0] = .5*(x[1]-x[0]); w[-1] = .5*(x[-1]-x[-2])
    w[1:-1] = .5*(x[2:]-x[:-2])
    return w


def xi_from_model(pf: ForestPower, cfg: Config, nk: int = 6400,
                  kmin: float = 1e-4, kmax: float = 35.0) -> XiTable:
    """Hankel/cosine transform of ``ForestPower`` on logarithmic k grids."""
    rp = np.arange(0, cfg.xi_max + .5*cfg.xi_step, cfg.xi_step)
    rz = rp.copy()
    kp = np.geomspace(kmin, kmax, nk)
    kz = np.geomspace(kmin, kmax, nk)
    wp = _trap_weights_log(kp) * kp**2/(2*np.pi)
    wz = _trap_weights_log(kz) * kz/np.pi  # 2 * int_0^inf dk/(2pi)
    J0 = j0(np.outer(rp, kp)); J1 = j1(np.outer(rp, kp))
    xi=np.zeros((len(rz),len(rp))); xirp=np.zeros_like(xi); variance=0.
    for i in range(0,nk,128):
        kzpart=kz[i:i+128]; wzpart=wz[i:i+128]
        P=pf(kzpart[:,None],kp[None,:])
        C=np.cos(np.outer(rz,kzpart))*wzpart
        xi+=C@((P*wp)@J0.T)
        xirp-=C@((P*(wp*kp))@J1.T)
        variance+=float(np.sum(P*wzpart[:,None]*wp))
    xi=xi.T; xirp=xirp.T
    return XiTable(rp, rz, xi, xirp,
                   {"provider": "model", "nk": nk, "variance_integral": variance,
                    "kmin": kmin, "kmax": kmax})


@njit(parallel=True,cache=True)
def _data_hist_kernel(pix_start, chi, delta, weight, a, b, theta,
                      rpmax, rzmax, step, slab, slab_edges, slab_index, region, chi_zedges):
    """Measured cells, split by pair type (index 0 = both pixels in region A, 1 = one in A and one in B) and by
    the pair's mean distance (bins ``chi_zedges``, the redshift bins of the evolution fit; pairs outside the
    edges are dropped). B x B pairs are dropped, because a region-B pixel also carries Lyb absorption from a much
    more distant slab and two of them correlate through it. Returns (num, den, chisum): the weighted products,
    the weights and the weighted mean-distance sum per cell, each [2, n_z, n_r, n_r]."""
    nr = int(np.ceil(rpmax/step)); nz = chi_zedges.size-1
    nums = np.zeros((24,2,nz,nr,nr),np.float64)
    dens = np.zeros((24,2,nz,nr,nr),np.float64)
    csum = np.zeros((24,2,nz,nr,nr),np.float64)
    for chunk in prange(24):
        for ip in range(chunk*a.size//24,(chunk+1)*a.size//24):
            aa, bb = a[ip], b[ip]
            q0=pix_start[bb]; qend=pix_start[bb+1]
            for p in range(pix_start[aa], pix_start[aa+1]):
                cp=float(chi[p])
                while q0<qend and float(chi[q0])<cp-rzmax: q0+=1
                q=q0
                while q<qend and float(chi[q])<=cp+rzmax:
                    rz = abs(float(chi[p])-float(chi[q]))
                    rp = .5*(float(chi[p])+float(chi[q]))*theta[ip]
                    cm=.5*(float(chi[p])+float(chi[q]))
                    nb = (1 if region[p]>0 else 0)+(1 if region[q]>0 else 0)
                    selected=slab[p]>=0 and slab[q]>=0 and nb<2
                    if slab_index>=0: selected=selected and slab_edges[slab_index,0]<=cm<slab_edges[slab_index,1]
                    if rz < rzmax and rp < rpmax and selected:
                        iz=-1
                        for k in range(nz):
                            if chi_zedges[k]<=cm<chi_zedges[k+1]:
                                iz=k; break
                        if iz>=0:
                            i, j = int(rp/step), int(rz/step)
                            ww = float(weight[p])*float(weight[q])
                            nums[chunk,nb,iz,i,j] += ww*float(delta[p])*float(delta[q])
                            dens[chunk,nb,iz,i,j] += ww
                            csum[chunk,nb,iz,i,j] += ww*cm
                    q+=1
    return nums.sum(axis=0),dens.sum(axis=0),csum.sum(axis=0)


def xi_from_data(sl: SightlineSet, cfg: Config) -> XiTable:
    """Measure, smooth and differentiate the observed cross-sightline xi."""
    try:
        from .pairs import find_pairs
    except ImportError:
        from pairs import find_pairs
    chi_min = max(1.0, float(np.min(sl.chi)))
    p = find_pairs(sl, cfg.xi_max/chi_min)
    region=np.asarray(getattr(sl,"region",None) if getattr(sl,"region",None) is not None
                      else np.zeros(len(sl.chi),np.int8),np.int8)
    from cosmo import chi as chi_of_z, z_of_chi
    # Redshift bins of the evolution fit (pairs are binned by mean distance). Without evolution, one bin that
    # spans every pair: the z edges are then only a bookkeeping device and no pair is dropped.
    zed = np.asarray(getattr(cfg, "xi_z_edges", None) or (0., 20.), float) if getattr(cfg, "xi_z_evolution", False) else np.array([0., 20.])
    chi_zedges = np.array([float(chi_of_z(z)) if z > 0 else 0. for z in zed]) if zed[-1] < 20 else np.array([0., 1e9])
    if zed[-1] < 20: chi_zedges[-1] = np.nextafter(chi_zedges[-1], np.inf)   # closed upper edge
    num, den, csum = _data_hist_kernel(sl.pix_start, sl.chi, sl.delta, sl.w,
                                 p[0], p[1], p[4], cfg.xi_max,
                                 cfg.xi_max, 1.0,sl.slab,
                                 np.asarray([[float(__import__("cosmo").chi(z)) for z in b] for b in cfg.slabs]),cfg.slab_index,
                                 region, chi_zedges)
    t = xi_from_counts(num.sum(axis=(0, 1)), den.sum(axis=(0, 1)), cfg, len(p[0]))
    # kept off `meta`, which is JSON-serialised on save
    t.counts_by_type = {"AA": (num[0].sum(axis=0), den[0].sum(axis=0)), "AB": (num[1].sum(axis=0), den[1].sum(axis=0))}
    # z-resolved cells for the evolution fit: [n_z, n_r, n_r] summed over pair types, and per type
    t.counts_z = {"all": (num.sum(axis=0), den.sum(axis=0), csum.sum(axis=0)),
                  "AA": (num[0], den[0], csum[0]), "AB": (num[1], den[1], csum[1]),
                  "chi_edges": chi_zedges, "z_edges": zed if zed[-1] < 20 else None}
    t.meta["pair_weight_by_type"] = {"AA": float(den[0].sum()), "AB": float(den[1].sum())}
    if zed[-1] < 20:
        dz = den.sum(axis=(0, 2, 3)); cz = csum.sum(axis=(0, 2, 3))
        t.meta["pair_weight_by_z"] = {"z_edges": zed.tolist(), "weight": dz.tolist(),
                                      "mean_z": [float(z_of_chi(c/w)) if w > 0 else None for c, w in zip(cz, dz)]}
    return t


def xi_from_counts(num,den,cfg,pair_count=0):
    """Smooth saved 1-Mpc/h counts, allowing reproducible development scans."""
    raw = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    # Normalised convolution prevents empty cells being interpreted as zeros.
    valid = (den > 0).astype(float)
    # Width is chosen on seeds 100-104, then frozen for every validation sample.
    smooth = gaussian_filter(raw*valid, cfg.xi_smoothing, mode="reflect")
    norm = gaussian_filter(valid, cfg.xi_smoothing, mode="reflect")
    coarse = np.divide(smooth, norm, out=np.zeros_like(smooth), where=norm > 1e-8)
    centers = np.arange(coarse.shape[0]) + .5
    fine = np.arange(0, cfg.xi_max + .5*cfg.xi_step, cfg.xi_step)
    from scipy.interpolate import RectBivariateSpline
    # Correlation is even in both separation coordinates. Extend bin centers
    # across zero before differentiating; clipping all r<0.5 to the first
    # center produced a spurious nonzero gradient at r_perp=0.
    signed=np.r_[-centers[::-1],centers]
    even=np.block([[coarse[::-1,::-1],coarse[::-1,:]],
                   [coarse[:,::-1],coarse]])
    spl = RectBivariateSpline(signed, signed, even, kx=3, ky=3, s=0)
    evaluation=np.minimum(fine,centers[-1])
    xi = spl(evaluation,evaluation)
    xirp = spl(evaluation,evaluation,dx=1)
    return XiTable(fine, fine, xi, xirp,
                   {"provider": "data", "coarse_step": 1.0,"smoothing_bins":cfg.xi_smoothing,"boundary":"even reflection about zero",
                    "accepted_weight": float(den.sum()), "n_sightline_pairs": pair_count},(num,den))


def direct_variance(pf, nk=700, kmin=1e-4, kmax=35.0):
    kp = np.geomspace(kmin, kmax, nk); kz = np.geomspace(kmin, kmax, nk)
    wp = _trap_weights_log(kp)*kp**2/(2*np.pi)
    wz = _trap_weights_log(kz)*kz/np.pi
    return float(np.sum(pf(kz[:,None], kp[None,:])*wz[:,None]*wp[None,:]))
