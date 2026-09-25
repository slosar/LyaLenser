"""Coordinate-shift bookkeeping tests (not a physical calibration)."""
from __future__ import annotations

import numpy as np
from lyalenser.config import SightlineSet
from lyalenser.pairs import find_pairs,accumulate,pair_midpoint_regions
from lyalenser.templates import Template,curl
from lyalenser.amplitude import amplitude


def shift_positions(sl,alpha,A):
    """Observed positions for an injection: theta' = theta - A alpha."""
    al=np.asarray(alpha,float); dec_rad=np.deg2rad(sl.dec)
    ra=sl.ra-np.rad2deg(A*al[:,0]/np.maximum(np.cos(dec_rad),1e-8))
    dec=sl.dec-np.rad2deg(A*al[:,1])
    return SightlineSet(sl.qid.copy(),ra,dec,sl.zq.copy(),sl.pix_start.copy(),
                        sl.chi.copy(),sl.delta.copy(),sl.w.copy(),sl.slab.copy(),
                        {**sl.attrs,"position_shift_amplitude":float(A),"position_shift_sign":"-alpha"},
                        sl.region.copy())


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


def paired_jackknife(x, values, samples, regions):
    """Paired slope uncertainty using matched delete-one regions across all shifts.

    Missing regions leave the corresponding full-sample estimate unchanged.
    The covariance between +A and -A is essential: individual amplitude errors
    cannot be propagated as though the shifted catalogues were independent.
    """
    union = np.unique(np.concatenate(regions))
    aligned = np.repeat(np.asarray(values, float)[:, None], len(union), axis=1)
    for i, (jk, reg) in enumerate(zip(samples, regions)):
        aligned[i, np.searchsorted(union, reg)] = jk
    slopes = np.array([paired_slopes(x, col)[0] for col in aligned.T])
    n = len(union)
    error = np.sqrt((n - 1) / n * np.sum((slopes - slopes.mean()) ** 2)) if n > 1 else np.nan
    return {'paired_slope_jk_error': float(error), 'slope_jk_samples': slopes.tolist(),
            'jk_regions': union.tolist()}


def injection_test(sl,xi_table,alpha_inj,A_list,cfg,templates=None,output=None,expectation=False):
    """Shift the sightline positions by -A alpha_inj, preserving A/B labels and production selection, and refit.

    ``expectation`` replaces the measured delta_p delta_q by the table's xi at the true (unshifted) separation
    (pairs.accumulate ``true_positions``): the noise-free expectation of the same test, including every boundary
    crossing of the r_perp cuts and the band-basis representation of the deflection. Unity is expected for a
    self-consistent lensable table and matching source-distance treatment; displacing an unlensed instrumental
    term while omitting its derivative from the response need not give a unit slope.
    """
    if templates is None:
        raise ValueError("injection_test requires science, curl and junk templates")
    vals=[]; curls=[]; errs=[]; samples=[]; regions=[]
    for A in A_list:
        shifted=shift_positions(sl,alpha_inj,A)
        ps=find_pairs(shifted,cfg.r_perp_max/max(float(shifted.chi.min()),1))
        cat=accumulate(shifted,ps,xi_table,cfg,true_positions=np.column_stack((sl.ra,sl.dec)) if expectation else None)
        reg=pair_midpoint_regions(cat,shifted,cfg.nside_jk)
        if cfg.nside_jk==8 and len(np.unique(reg))<30:
            reg=pair_midpoint_regions(cat,shifted,16)
        r=amplitude(cat,templates,cfg.g1,reg)
        from lyalenser.amplitude import common_science
        from lyalenser.amplitude import curl_amplitude
        science=common_science(r)
        vals.append(science["A"]); curls.append(curl_amplitude(r)[0])
        samples.append(science['jk']); regions.append(r.regions)
        errs.append(r.jk_error.tolist())
        if output is not None:
            tag="expectation" if expectation else "injection"
            r.save(output,f"{tag}/fit_{A}")
            cat.save(output,f"{tag}/catalogue_{A}")
    x=np.asarray(A_list,float); y=np.asarray(vals)
    slope,per=paired_slopes(x,y)
    cslope=np.polyfit(x,np.asarray(curls),1)[0] if len(x)>1 else np.nan
    return {"A_injected":x,"A_hat":y,"curl":np.asarray(curls),"errors":np.asarray(errs),
            **paired_jackknife(x,vals,samples,regions),
            "paired_slope":float(slope),"paired_slopes_by_amplitude":per,"curl_slope":float(cslope),
            "expectation":bool(expectation),
            "description":"coordinate bookkeeping test (noise-free expectation)" if expectation else "coordinate bookkeeping test; not a physical calibration"}
