"""The fitted forest correlation table of a production run: reading a saved table, and the fit of the
measured 1 Mpc/h pair counts to the projected model basis with the spline correction (`xi_spline`).

The response kernel of the pair estimator is G = chi d xi / d r_perp of this table; see docs/pipeline.md.
"""
import json
import numpy as np
from lyalenser.xi_model import XiTable, xi_from_data


def read_xi(path, group='xi'):
    """Read an `XiTable` saved by `XiTable.save` (layered tables carry `chi_nodes`)."""
    import h5py
    with h5py.File(path) as f:
        g = f[group]
        return XiTable(*(g[k][()] for k in ('r_perp', 'r_par', 'xi', 'xi_rp')), json.loads(g.attrs['meta']),
                       chi_nodes=(g['chi_nodes'][()] if 'chi_nodes' in g else None))


def correction_for(basis, cfg):
    """The production spline correction (`xi_spline.XiCorrection`) on the envelope of the basis' reference shape.
    `cfg.xi_knots` selects the knot set: 'bicubic' (16 + 6 coefficients, fiducial), 'medium' (30 + 6) or
    'fine' (42 + 6); `cfg.xi_same_wavelength` switches the same-wavelength term."""
    from lyalenser.xi_spline import Envelope, XiCorrection
    from lyalenser.xi_fit import BASIS
    b = basis['projected']; c = [1., 2 * 1.4, 1.4 ** 2]
    ref = XiTable(b[BASIS[0]].r_perp, b[BASIS[0]].r_par, sum(x * b[k].xi for x, k in zip(c, BASIS)),
                  sum(x * b[k].xi_rp for x, k in zip(c, BASIS)), {})
    KN = {'bicubic': ((3., 30.), (0., 30.)),
          'medium': ((3., 8., 16., 30.), (0., 6., 30.)),
          'fine': ((3., 6., 10., 16., 30.), (0., 4., 10., 30.))}[getattr(cfg, 'xi_knots', 'bicubic')]
    return XiCorrection(Envelope(ref), rp_knots=KN[0], rz_knots=KN[1],
                        rz_sw=1. if getattr(cfg, 'xi_same_wavelength', True) else 0.)


def table_for(sl, cfg, basis, counts=None):
    """Measured 1 Mpc/h counts of a sightline set -> fitted response table (`xi_fit.FittedTable`).

    `cfg.xi_correction` selects the two-parameter Kaiser fit of the projected model basis ('none') or the same
    fit plus the spline correction and the same-wavelength term ('spline', the production choice)."""
    num, den = (xi_from_data(sl, cfg).counts if counts is None else counts)
    if getattr(cfg, 'xi_correction', 'none') == 'spline':
        from lyalenser.xi_spline import fit_corrected_table
        from lyalenser.xi_fit import FittedTable
        corr = correction_for(basis, cfg)
        tab, par, _ = fit_corrected_table(num, den, basis['projected'], basis['coarse'], cfg, corr, ridge=cfg.xi_correction_ridge)
        return FittedTable(tab, par)
    from lyalenser.xi_fit import fit_model_table
    return fit_model_table(num, den, basis['projected'], cfg, basis['coarse'])
