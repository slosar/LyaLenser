"""Joint response fit over all slice templates (iteration 14, user request ii).

The production fits treat one template at a time (an 11 x 11 response matrix per slice: five science bands, five
curl partners, the junk band) and the combined template separately. Here ALL slice templates enter one response
matrix (5 slices x 11 = 55 components for the fiducial bands): q_i = sum_pairs G_i (dd - xi) xi', F_ij =
sum_pairs G_i G_j xi'^2 (with the source-distance terms of `amplitude._partials`, from the templates' derivative maps), accumulated PER JACKKNIFE REGION
without the (N_t, N_t, N_pair) array (7 MB instead of 124 GB). From (q, F, mf) and their region partials:
  * the global common science amplitude (all slice-band science components share one amplitude; every curl and
    junk component is a free nuisance),
  * the per-slice common amplitudes fitted JOINTLY (one amplitude per slice, others marginalised),
  * the per-band common amplitudes fitted jointly (one amplitude per band across slices),
each with the jackknife error and samples, and the normalised response matrix for the figure.
"""
from __future__ import annotations
import numpy as np

SCIENCE_KINDS = {"signal", "truth", "injection", "response", "random"}


def partials_by_region(cat, templates, g1, regions, bins=None):
    """(names, region ids, pq [nr, nt], pF [nr, nt, nt], pmf [nr, nt]) with per-region matrix products."""
    # the same contraction as amplitude._partials (d, dp, sp per template: alpha difference, derivative-map
    # difference and sum across the pair), accumulated per region with matrix products. The pair scalars are
    # formed region by region: for 55 components and 1e7 pairs the full (N_t, N_pair) arrays would need ~17 GB.
    from lyalenser.amplitude import _alpha_name, _derivative
    names = [_alpha_name(t, i)[1] for i, t in enumerate(templates)]; nt = len(names)
    A = np.stack([_alpha_name(t, i)[0] for i, t in enumerate(templates)])          # (nt, Nq, 2)
    Ad = np.stack([_derivative(t, g1) for t in templates])                             # (nt, Nq, 2)
    bins = np.arange(6) if bins is None else np.atleast_1d(bins)
    x = cat.accum[:, :, bins].sum(axis=2)
    regvals = np.unique(regions); nr = len(regvals); inv = np.searchsorted(regvals, regions)
    pq = np.zeros((nr, nt)); pmf = np.zeros((nr, nt)); pF = np.zeros((nr, nt, nt))
    order = np.argsort(inv, kind='stable'); bounds = np.searchsorted(inv[order], np.arange(nr + 1))
    thx = cat.thx.astype(float); thy = cat.thy.astype(float)
    for k in range(nr):
        idx = order[bounds[k]:bounds[k + 1]]
        if len(idx) == 0: continue
        ia = cat.a[idx]; ib = cat.b[idx]; tx = thx[idx]; ty = thy[idx]
        da = A[:, ia] - A[:, ib]; dda = Ad[:, ia] - Ad[:, ib]; dsa = Ad[:, ia] + Ad[:, ib]
        D = tx * da[:, :, 0] + ty * da[:, :, 1]; Dp = tx * dda[:, :, 0] + ty * dda[:, :, 1]; Sp = tx * dsa[:, :, 0] + ty * dsa[:, :, 1]
        xk = x[idx]
        pq[k] = D @ xk[:, 0] + Dp @ xk[:, 1] + Sp @ xk[:, 2]; pmf[k] = D @ xk[:, 8] + Dp @ xk[:, 9] + Sp @ xk[:, 10]
        DDp = D * xk[:, 4]; DS = D * (.5 * xk[:, 6])
        pF[k] = (D * xk[:, 3]) @ D.T + (DDp @ Dp.T + Dp @ DDp.T) + (Dp * xk[:, 5]) @ Dp.T + (DS @ Sp.T + Sp @ DS.T) + (Sp * xk[:, 7]) @ Sp.T
    return names, regvals, pq, pF, pmf


class JointFit:
    def __init__(self, names, kinds, groups, regvals, pq, pF, pmf):
        """``groups``: per component, the group label of its science amplitude (None for nuisance)."""
        self.names = list(names); self.kinds = list(kinds); self.groups = list(groups); self.regvals = np.asarray(regvals)
        self.pq = pq; self.pF = pF; self.pmf = pmf; self.q = pq.sum(axis=0); self.F = pF.sum(axis=0); self.mf = pmf.sum(axis=0)

    def _design(self, science_groups):
        """Columns: one per science group (its components share the amplitude), then one per nuisance component."""
        nt = len(self.names); cols = []
        for g in science_groups:
            c = np.zeros(nt); c[[i for i, gg in enumerate(self.groups) if gg == g]] = 1.; cols.append(c)
        for i, gg in enumerate(self.groups):
            if gg is None: c = np.zeros(nt); c[i] = 1.; cols.append(c)
        return np.column_stack(cols)

    def fit(self, science_groups):
        """Amplitudes of the science groups with every nuisance free; jackknife over the regions."""
        M = self._design(science_groups); ng = len(science_groups)
        def solve(F, y):
            Fr = M.T @ F @ M; return np.linalg.solve(Fr, M.T @ y)
        A = solve(self.F, self.q - self.mf)[:ng]; Finv = np.linalg.inv(M.T @ self.F @ M); sig = np.sqrt(np.diag(Finv))[:ng]
        nr = len(self.regvals); jk = np.zeros((nr, ng))
        for r in range(nr):
            jk[r] = solve(self.F - self.pF[r], (self.q - self.pq[r]) - (self.mf - self.pmf[r]))[:ng]
        dif = jk - jk.mean(axis=0); cov = (nr - 1) / nr * dif.T @ dif
        return {'groups': list(science_groups), 'A': A.tolist(), 'sigma_F': sig.tolist(), 'jk_error': np.sqrt(np.diag(cov)).tolist(), 'jk_cov': cov.tolist(), 'jk_samples': jk, 'regions': self.regvals}

    def normalised_response(self):
        d = np.sqrt(np.diag(self.F)); return self.F / np.outer(d, d)


def build_joint(cat, slice_templates: dict, g1, regions, bins=None):
    """``slice_templates``: {slice name: list of Template objects (science, curl, junk)}. Returns a JointFit whose
    component groups are (slice, band) for science components and None for curl/junk."""
    templates = []; kinds = []; groups = []; labels = []
    for sname, tl in slice_templates.items():
        for t in tl:
            templates.append(t); kinds.append(t.kind); labels.append(f'{sname}:{t.name}')
            groups.append((sname, t.name) if t.kind in SCIENCE_KINDS else None)
    names, regvals, pq, pF, pmf = partials_by_region(cat, templates, g1, regions, bins)
    return JointFit(labels, kinds, groups, regvals, pq, pF, pmf)


def block_diagonal(jf: JointFit):
    """The same fit with the response matrix restricted to its per-slice blocks (the 'diagonal R' variant)."""
    sl = np.array([n.split(':')[0] for n in jf.names]); keep = (sl[:, None] == sl[None, :]).astype(float)
    return JointFit(jf.names, jf.kinds, jf.groups, jf.regvals, jf.pq, jf.pF * keep[None, :, :], jf.pmf)


def drop_kinds(jf: JointFit, kinds=('curl',)):
    """The same fit with the components of the given kinds REMOVED from the model (their amplitudes fixed to
    zero instead of marginalised): the 'no curl' variant of the fiducial (user request, 2026-09-21)."""
    keep = np.array([k not in kinds for k in jf.kinds])
    return JointFit([n for n, k in zip(jf.names, keep) if k], [x for x, k in zip(jf.kinds, keep) if k], [g for g, k in zip(jf.groups, keep) if k],
                    jf.regvals, jf.pq[:, keep], jf.pF[:, keep][:, :, keep], jf.pmf[:, keep])


def standard_fits(jf: JointFit):
    """The global amplitude, the per-slice and the per-band amplitudes, all joint."""
    sci = [g for g in jf.groups if g is not None]; slices = list(dict.fromkeys(g[0] for g in sci)); bands = list(dict.fromkeys(g[1] for g in sci))
    out = {}
    j2 = JointFit(jf.names, jf.kinds, ['all' if g is not None else None for g in jf.groups], jf.regvals, jf.pq, jf.pF, jf.pmf); out['global'] = j2.fit(['all'])
    j2 = JointFit(jf.names, jf.kinds, [g[0] if g is not None else None for g in jf.groups], jf.regvals, jf.pq, jf.pF, jf.pmf); out['per_slice'] = j2.fit(slices)
    j2 = JointFit(jf.names, jf.kinds, [g[1] if g is not None else None for g in jf.groups], jf.regvals, jf.pq, jf.pF, jf.pmf); out['per_band'] = j2.fit(bands)
    # the curl null: all curl components share one amplitude, science components per (slice, band) free
    groups_curl = [g if g is not None else (('curl',) if k == 'curl' else None) for g, k in zip(jf.groups, jf.kinds)]
    j2 = JointFit(jf.names, jf.kinds, groups_curl, jf.regvals, jf.pq, jf.pF, jf.pmf); out['curl'] = j2.fit([('curl',)] + sci)
    out['curl'] = {k: (v[:1] if isinstance(v, list) and k in ('A', 'sigma_F', 'jk_error') else v) for k, v in out['curl'].items()}
    if 'curl' in jf.kinds:
        nc = drop_kinds(jf, ('curl',)); j2 = JointFit(nc.names, nc.kinds, ['all' if g is not None else None for g in nc.groups], nc.regvals, nc.pq, nc.pF, nc.pmf)
        out['no_curl_global'] = j2.fit(['all'])
        j2 = JointFit(nc.names, nc.kinds, [g[0] if g is not None else None for g in nc.groups], nc.regvals, nc.pq, nc.pF, nc.pmf); out['no_curl_per_slice'] = j2.fit(slices)
        j2 = JointFit(nc.names, nc.kinds, [g[1] if g is not None else None for g in nc.groups], nc.regvals, nc.pq, nc.pF, nc.pmf); out['no_curl_per_band'] = j2.fit(bands)
        nn = drop_kinds(jf, ('curl', 'junk')); j2 = JointFit(nn.names, nn.kinds, ['all' if g is not None else None for g in nn.groups], nn.regvals, nn.pq, nn.pF, nn.pmf)
        out['no_curl_no_junk_global'] = j2.fit(['all'])
    return out
