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


def paired_slopes(x,y):
    """Odd component per injected |A| ((y(+a)-y(-a))/2a) and the amplitude-weighted combination."""
    x=np.asarray(x,float); y=np.asarray(y,float); per={}
    num=den=0.0
    for a in sorted(set(abs(x[x!=0]))):
        if np.any(x==a) and np.any(x==-a):
            d=(y[x==a][0]-y[x==-a][0])/(2*a); per[float(a)]=float(d)
            num += a*(y[x==a][0]-y[x==-a][0]); den += 2*a*a
    slope=num/den if den else np.polyfit(x,y,1)[0]
    return float(slope),per


def injection_test(sl,xi_table,alpha_inj,A_list,cfg,templates=None,output=None,expectation=False):
    """Shift the sightline positions by -A alpha_inj, rebuild the pairs with the production selection and refit.

    ``expectation`` replaces the measured delta_p delta_q by the table's xi at the true (unshifted) separation
    (pairs.accumulate ``true_positions``): the noise-free expectation of the same test, including every boundary
    crossing of the r_perp cuts and the band-basis representation of the exact deflection. Its odd slope must
    approach 1 as A -> 0 if the bookkeeping (signs, metric, cuts, templates) is right (GATES v6).
    """
    if templates is None:
        raise ValueError("injection_test requires the full seven-component map basis")
    vals=[]; curls=[]; errs=[]
    for A in A_list:
        shifted=shift_positions(sl,alpha_inj,A)
        ps=find_pairs(shifted,cfg.r_perp_max/max(float(shifted.chi.min()),1))
        cat=accumulate(shifted,ps,xi_table,cfg,true_positions=np.column_stack((sl.ra,sl.dec)) if expectation else None)
        reg=pair_midpoint_regions(cat,shifted,cfg.nside_jk)
        if cfg.nside_jk==8 and len(np.unique(reg))<30:
            reg=pair_midpoint_regions(cat,shifted,16)
        r=amplitude(cat,templates,cfg.g1,reg)
        from run_mock_validation import common_science
        from amplitude import curl_amplitude
        vals.append(common_science(r)["A"]); curls.append(curl_amplitude(r)[0])
        errs.append(r.jk_error.tolist())
        if output is not None:
            tag="expectation" if expectation else "injection"
            r.save(output,f"{tag}/fit_{A}")
            cat.save(output,f"{tag}/catalogue_{A}")
    x=np.asarray(A_list,float); y=np.asarray(vals)
    slope,per=paired_slopes(x,y)
    cslope=np.polyfit(x,np.asarray(curls),1)[0] if len(x)>1 else np.nan
    return {"A_injected":x,"A_hat":y,"curl":np.asarray(curls),"errors":np.asarray(errs),
            "paired_slope":float(slope),"paired_slopes_by_amplitude":per,"curl_slope":float(cslope),
            "expectation":bool(expectation),
            "description":"coordinate bookkeeping test (noise-free expectation)" if expectation else "coordinate bookkeeping test; not a physical calibration"}
