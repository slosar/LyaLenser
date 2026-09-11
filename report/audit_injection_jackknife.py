"""Prewritten audit of the fitted curl slope; no mock generation or tuning."""
from pathlib import Path
import json,hashlib
import h5py
import numpy as np

root=Path('/data/LyaLenser/mocks/iteration3')
registration=json.loads((root/'injection_audit_registration.json').read_text())
assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==registration['code_sha256']
x=np.array([-2.,-1.,0.,1.,2.]); xc=x-x.mean()
fits=[]
with h5py.File(root/'injection.h5','r') as f:
    for a in x:
        g=f[f'injection/fit_{int(a)}']
        fits.append({key:g[key][:] for key in ('A','regions','jk_samples','jk_cov')})
regions=sorted(set().union(*(set(s['regions'].tolist()) for s in fits)))
full=np.array([s['A'][3:6].mean() for s in fits])
slope=float(xc@full/(xc@xc))
leave=[]
for region in regions:
    values=[]
    for fit in fits:
        where=np.flatnonzero(fit['regions']==region)
        # One shared HEALPix partition; an empty region leaves a fit unchanged.
        pars=fit['jk_samples'][where[0]] if len(where) else fit['A']
        values.append(pars[3:6].mean())
    leave.append(float(xc@np.array(values)/(xc@xc)))
leave=np.array(leave); n=len(leave)
error=float(np.sqrt((n-1)/n*np.sum((leave-leave.mean())**2)))
rms=float(np.sqrt(np.mean([np.diag(s['jk_cov'])[3:6] for s in fits])))
result={'method':'OLS curl slope on absolute mean-curl fits; common union of saved midpoint-HEALPix regions; standard delete-one-region jackknife',
        'registration':registration,
        'A_injected':x.tolist(),'absolute_mean_curl':full.tolist(),'curl_slope':slope,
        'joint_slope_jk_error':error,'frozen_rms_individual_band_jk_scale':rms,
        'regions':regions,'leave_region_slopes':leave.tolist(),
        'within_joint_one_sigma':bool(np.isfinite(error) and abs(slope)<=error),
        'within_frozen_rms_scale':bool(np.isfinite(rms) and abs(slope)<=rms),
        'note':'Additional uncertainty audit; the frozen acceptance implementation and its decision are retained separately.'}
(root/'injection_joint_jackknife_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('regions','leave_region_slopes')},indent=2))
