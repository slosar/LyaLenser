"""Coordinate-shift bookkeeping tests (not a physical calibration)."""
from __future__ import annotations

import numpy as np
try:
    from .config import SightlineSet
    from .pairs import find_pairs,accumulate,pair_midpoint_regions
    from .templates import Template,curl
    from .amplitude import amplitude
except ImportError:
    from config import SightlineSet
    from pairs import find_pairs,accumulate,pair_midpoint_regions
    from templates import Template,curl
    from amplitude import amplitude


def shift_positions(sl,alpha,A):
    """Observed positions for an injection: theta' = theta - A alpha."""
    al=np.asarray(alpha,float); dec_rad=np.deg2rad(sl.dec)
    ra=sl.ra-np.rad2deg(A*al[:,0]/np.maximum(np.cos(dec_rad),1e-8))
    dec=sl.dec-np.rad2deg(A*al[:,1])
    return SightlineSet(sl.qid.copy(),ra,dec,sl.zq.copy(),sl.pix_start.copy(),
                        sl.chi.copy(),sl.delta.copy(),sl.w.copy(),sl.slab.copy(),
                        {**sl.attrs,"position_shift_amplitude":float(A),"position_shift_sign":"-alpha"})


def injection_test(sl,xi_table,alpha_inj,A_list,cfg,junk_alpha=None):
    vals=[]; curls=[]; errs=[]
    t=Template(alpha_inj,"injection","injection")
    tc=Template(curl(alpha_inj),"injection_curl","curl")
    if junk_alpha is None:
        # Deterministic outside-template nuisance direction for bookkeeping
        # callers that do not own a map. Production validation passes the
        # actual Fourier junk band explicitly.
        phase=2*np.pi*np.arange(len(alpha_inj))/max(len(alpha_inj),1)
        scale=max(float(np.std(alpha_inj)),1e-8)
        junk_alpha=scale*np.column_stack((np.cos(phase),np.sin(3*phase)))
    tj=Template(junk_alpha,"junk","junk")
    for A in A_list:
        shifted=shift_positions(sl,alpha_inj,A)
        ps=find_pairs(shifted,cfg.r_perp_max/max(float(shifted.chi.min()),1))
        cat=accumulate(shifted,ps,xi_table,cfg)
        reg=pair_midpoint_regions(cat,shifted,cfg.nside_jk)
        r=amplitude(cat,[t,tc,tj],cfg.g1,reg)
        vals.append(r.A[0]); curls.append(r.A[1]); errs.append(r.jk_error.tolist())
    x=np.asarray(A_list,float); y=np.asarray(vals)
    # Paired +/- differences are exactly the odd component.
    num=den=0.0
    for a in sorted(set(abs(x[x!=0]))):
        if np.any(x==a) and np.any(x==-a):
            num += a*(y[x==a][0]-y[x==-a][0]); den += 2*a*a
    slope=num/den if den else np.polyfit(x,y,1)[0]
    cslope=np.polyfit(x,np.asarray(curls),1)[0] if len(x)>1 else np.nan
    return {"A_injected":x,"A_hat":y,"curl":np.asarray(curls),"errors":np.asarray(errs),
            "paired_slope":float(slope),"curl_slope":float(cslope),
            "description":"coordinate bookkeeping test; not a physical calibration"}
