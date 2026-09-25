"""Summarise the reproducible 100-draw null check and paper numerical corrections."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from lyalenser.cosmo import chi


def statistics(draws):
    A=np.array([d['A'] for d in draws]); err=np.array([d['jk_error'] for d in draws])
    n=len(A); scatter=A.std(ddof=1); sem=scatter/np.sqrt(n); rms=np.sqrt(np.mean(err**2))
    rng=np.random.default_rng(9173); idx=rng.integers(n,size=(10000,n))
    ratios=A[idx].std(axis=1,ddof=1)/np.sqrt(np.mean(err[idx]**2,axis=1))
    return {'n':n,'mean':float(A.mean()),'sem':float(sem),'scatter':float(scatter),
            'rms_jk_error':float(rms),'mean_t':float(A.mean()/sem),
            'mean_p_two_sided':float(2*student_t.sf(abs(A.mean()/sem),n-1)),
            'scatter_over_rms_jk':float(scatter/rms),
            'scatter_over_rms_jk_bootstrap95':np.quantile(ratios,[.025,.975]).tolist(),
            'normalised_mean':float(np.mean(A/err)),
            'normalised_scatter':float(np.std(A/err,ddof=1))}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--directory',type=Path,required=True)
    ap.add_argument('--out',type=Path,default=ROOT/'results/random_template_null_100_v4.json')
    args=ap.parse_args()
    out={'seed':2026,'fit':'11-component combined-template diagnostic, not the 55-component fiducial slice fit',
         'bootstrap_draws':10000,'statistics':{}}
    for stat,old in [('auto','dr1_lowz_v7d'),('cross','dr1_qso_v1d')]:
        draws=[]
        for path in sorted(args.directory.glob(f'null_{stat}*.json')):
            shard=json.loads(path.read_text())
            assert shard['seed']==2026 and shard['statistic']==stat
            draws.extend(shard['draws'])
        draws.sort(key=lambda d:d['index'])
        if [d['index'] for d in draws]!=list(range(100)): raise ValueError('missing/duplicate draws')
        oldA=np.array(json.loads((ROOT/'results'/f'{old}.json').read_text())['random_templates']['amplitudes'])
        difference=np.max(np.abs(np.array([d['A'] for d in draws[:20]])-oldA))
        np.testing.assert_allclose([d['A'] for d in draws[:20]],oldA,rtol=2e-5,atol=2e-6)
        rec={'all100':statistics(draws),'original20':statistics(draws[:20]),
             'additional80':statistics(draws[20:]),'max_original20_amplitude_difference':float(difference),
             'draws':draws}
        out['statistics'][stat]=rec
        print(stat,json.dumps({k:v for k,v in rec.items() if k!='draws'},indent=1))
    auto=json.loads((ROOT/'results/dr1_lowz_v7d.json').read_text())
    z=np.linspace(1.96,3.,1001); zb=(1+z)*1215.67/1025.72-1; separation=chi(zb)-chi(z)
    out['paper_numerical_checks']={'sightline_density_deg2':auto['forests']/auto['area_deg2_nside64'],
                                  'lya_lyb_separation_mpch':[float(separation.min()),float(separation.max())]}
    args.out.write_text(json.dumps(out,indent=1)+'\n')


if __name__=='__main__': main()
