"""Format saved iteration-3 results; never run mocks or change acceptance decisions.

Run after code/pipeline/run_mock_validation.py --phase rebuild.
"""
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def number(x):
    return f"{x:.6g}"


def stat(s, unit="A"):
    return (f"{number(s['mean'])} ± {number(s['sem'])} SEM; "
            f"95% residual bound {number(s['bound95'])} {unit}; N={s['n']}")


def cell(text):
    return str(text).replace("|", "&#124;").replace("\n", " ")


def measurement(row):
    s = row['measured']
    name = row['acceptance']
    if 'mean' in s:
        result = stat(s, "(slope units)" if 'slope' in name else "A")
        if 'prediction' in s:
            result += "; predicted " + stat(s['prediction'])
            result += "; observed − predicted " + stat(s['prediction_difference'])
        if 'raw' in s:
            result += "; raw " + stat(s['raw'])
        if 'scatter_over_rms_jk' in s:
            result += f"; scatter/RMS JK={number(s['scatter_over_rms_jk'])}; scatter/RMS σF={number(s['scatter_over_rms_sigma_F'])}"
        return result
    if 'fixed_absolute_slope' in s:
        d = s['diagnostic_difference']
        return ("Fixed slope " + stat(s['fixed_absolute_slope'], "(slope units)")
                + "; refitted slope " + stat(s['refitted_absolute_slope'], "(slope units)")
                + "; diagnostic refitted − fixed slope " + stat(d, "(slope units)")
                + "; diagnostic A_true=1 difference " + stat(d['absolute_A1_refit_minus_fixed']))
    if 'p_value' in s:
        return f"Hotelling T²={number(s['T2'])}; F={number(s['F'])}; df={s['df']}; p={number(s['p_value'])}"
    if 'ratios' in s:
        return "; ".join(f"{n}: {number(v)} ± {number(e)} SEM" for n, v, e in zip(('klkl','klkc','kckc'), s['ratios'], s['sem']))
    if 'benchmark' in s:
        b = s['benchmark']
        return f"{number(b['pixel_pairs_per_s'])} pixel pairs/s; {b['threads']} threads; peak {number(s['peak_gib'])} GiB ({number(s['peak_gib']*1024**3/1e9)} GB)"
    if 'passed' in s:
        return f"{s['passed']} passed, {s['failed']} failed; {number(s['wall_s'])} s"
    if 'regression_stderr' in s:
        return f"Slope {number(s['slope'])}; descriptive regression error {number(s['regression_stderr'])}; ensemble SEM and 95% bound unavailable (one realization); absolute A={s['A']}"
    if 'flags' in s:
        return "; ".join(f"{k} off: absolute A={number(v['absolute_A'])}, diagnostic shift={number(v['diagnostic_difference'])}" for k,v in s['flags'].items()) + "; nested margins persisted"
    if name == 'injection bookkeeping':
        return f"Science odd slope {number(s['slope'])}; curl slope {number(s['curl_slope'])}; frozen RMS individual-band JK error scale {number(s['curl_error'])}; joint slope-error audit below"
    return json.dumps(s, ensure_ascii=False)


def main():
    result = json.loads((HERE/'mock_validation.json').read_text())
    assert result['iteration'] == 3
    rows = result['acceptance']
    assert len(rows) == 30, 'The complete prospective table must contain all 30 gates.'
    seeds = result['seeds']
    rec = [s['seed'] for s in seeds if s['role'] in ('recovery','extension')]
    null = [s['seed'] for s in seeds if s['role'] in ('null','extension')]
    failures = [r for r in rows if not r['pass']]
    lines = ['# Stage A iteration 3 validation', '',
             f"**{len(rows)-len(failures)}/{len(rows)} gates PASS. " + ('Stage B remains blocked.' if failures else 'All declared gates pass.') + '**', '',
             'Recovery and null amplitudes are absolute `F^-1(q-mf)`. Each mean-field subtraction belongs to the sample being fitted. Paired changes are labelled diagnostics.', '',
             'The scalar common-science estimate constrains the three science-band amplitudes to one value while fitting the three curls and junk as separate nuisance components. The full seven-component response and amplitude vector are also saved.', '',
             'The reported bound is `abs(mean − target) + t(0.975, N−1) × SEM`, with sample SEM. Amplitude bounds have A units; slope bounds are dimensionless. Normalization decisions use the predeclared ±0.03 mean-slope tolerance, with uncertainty also reported. A two-SEM consistency alone does not establish precise normalization.', '',
             f'Recovery seeds ({len(rec)}): {rec}.', '', f'Null seeds ({len(null)}): {null}.', '',
             '| Gate | Measurement | Frozen tolerance | Result |', '|---|---|---|---|']
    for row in rows:
        lines.append('| ' + ' | '.join(map(cell, (row['acceptance'], measurement(row), row['tolerance'], 'PASS' if row['pass'] else 'FAIL'))) + ' |')
    root = Path(seeds[0]['files']['mock']).parent
    decision = json.loads((root/'extension_decision.json').read_text())
    result['supplementary_audits'] = {
        name: json.loads((root/name).read_text())
        for name in ('final_persistence_audit.json', 'control_persistence_audit.json', 'random_stream_audit.json', 'final_diagnostic_summary.json', 'injection_joint_jackknife_audit.json')
        if (root/name).exists()
    }
    (HERE/'mock_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    lines += ['', '## Prospective stopping decision', '',
              'After the initial 40 seeds, recovery: ' + stat(decision['recovery']) + '; response-only null: ' + stat(decision['null']) + '.', '',
              f"Extension triggered: {decision['extend']}. Added seeds: {decision['seeds']}. The decision was made once; no seeds beyond 59 were used.", '',
              '## Absolute amplitude points and shape', '',
              'All per-seed amplitude-grid points, raw/corrected fits, jackknife errors, predicted responses, and shape estimates are retained in [the JSON report](mock_validation.json). Figures show ensemble SEM, except individual combined-recovery points, which show midpoint-jackknife errors.', '',
              '- [Physical normalization](figures/mock_physical_normalisation.pdf)',
              '- [Absolute combined recovery](figures/mock_combined_ensemble.pdf)',
              '- [Six-bin shape](figures/mock_shape.pdf)',
              '- [Spectra](figures/mock_spectra.pdf)', '',
              'Shape bins 0–5 are, in order, transverse separations 0–10, 10–20, and 20–30 Mpc/h, each split into radial separations 0–10 and 10–30 Mpc/h.', '',
              '## Same-realization margin diagnostic', '',
              'All margins below use seed 0 with the same larger density realization, forest, sightlines, and CMB map. Differences from the 150 Mpc/h selection are diagnostics. Flag and margin checks each use one realization, so an ensemble SEM or 95% bound is unavailable.', '',
              '| Margin (Mpc/h) | Absolute A | Diagnostic change from 150 |', '|---|---|---|']
    margins = next(s for s in seeds if s['seed'] == 0)['margin']
    for margin, value in sorted(margins.items(), key=lambda kv: float(kv[0])):
        lines.append(f"| {number(float(margin))} | {number(value)} | {number(value-margins['150.0'])} |")
    lines += ['', '## Coordinate injection bookkeeping', '',
              'The values below are absolute fits on one realization. On this symmetric injection grid, the reported odd slope equals the least-squares slope with a free intercept. This tests coordinate bookkeeping; physical normalization is assessed by the fresh remapping ensemble. An ensemble SEM or 95% bound is unavailable for this single-realization control.', '',
              '| Injected A | Absolute fitted A | Absolute mean curl amplitude |', '|---|---|---|']
    injection = result['extras']['injection']
    for x, y, curl in zip(injection['A_injected'], injection['A_hat'], injection['curl']):
        lines.append(f"| {number(x)} | {number(y)} | {number(curl)} |")
    joint = result['supplementary_audits'].get('injection_joint_jackknife_audit.json')
    if joint:
        lines += ['', f"Joint curl-slope jackknife audit: slope {number(joint['curl_slope'])}, slope error {number(joint['joint_slope_jk_error'])}; within one sigma: {'PASS' if joint['within_joint_one_sigma'] else 'FAIL'}. The frozen RMS individual-band error scale is {number(joint['frozen_rms_individual_band_jk_scale'])}. The audit uses the common union of saved midpoint-HEALPix regions and applies the original jackknife variance formula to the leave-region slopes. It is reported separately from the frozen acceptance decision."]
    lines += ['',
              '## Wall time and memory', '', '```json', json.dumps(result['timing'], indent=2), '```', '',
              f"Successful campaign peak RSS: {number(result['peak_memory_gb'])} GiB ({number(result['peak_memory_gb']*1024**3/1e9)} GB).", '',
              'The initial development process was interrupted with SIGKILL; its last observed RSS is not its known termination peak. The allocation repair preceded validation. Overlapping development-job durations must not be summed as elapsed campaign time.', '',
              '## Provenance and reconstruction', '',
              'Numerical choices were frozen on development seeds 100–104. [GATES.md](../GATES.md) predates validation; its hash, source hashes, configuration hash and library versions are recorded below. The reporting formatter does not recompute fits or change PASS/FAIL decisions.', '',
              '```json', json.dumps(result['provenance'], indent=2), '```', '',
              'From `code/pipeline/`, run `python run_mock_validation.py --phase rebuild`; then, from the repository root, run `python report/audit_injection_jackknife.py` and `python report/render_iteration3.py`. This uses the persisted ensemble, controls and diagnostics without generating new mocks.', '',
              '## Remaining failures', '']
    for row in failures:
        lines += [f"- **{row['acceptance']}**: {row['detail'] or 'The frozen tolerance was not met.'}"]
    if not failures:
        lines.append('None.')
    if 'random_stream_audit.json' in result['supplementary_audits']:
        lines += ['', '### Additional unresolved independence issue', '',
                  'A source audit after the freeze found that the field and CMB-noise generators initialize separate random generators with the same seed. Their first normal draws reuse the same stream prefix. The saved small-case reproduction confirms exact prefix reuse. The effect after Fourier filtering and array reshaping has not been isolated; it is not claimed to explain the measured offsets.', '',
                  'This issue is outside the 17 listed findings and is not a post-hoc change to the acceptance table. The scientific source and ensemble were kept frozen. A new development/freeze/validation cycle with independent random streams is required before the stochastic results can certify independent instrumental noise. The audit is included in the JSON report.']
    lines += ['', 'Numerical diagnoses and implementation limitations are recorded in the iteration-3 section of [NOTES.md](../code/pipeline/NOTES.md). No failed gate was repaired by changing validation seeds or tolerances.', '']
    (HERE/'mock_validation.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
