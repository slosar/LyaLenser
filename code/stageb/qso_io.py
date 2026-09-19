"""DR1 quasar catalogue for the quasar x forest cross-correlation (iteration 12).

`QSO_cat_iron_cumulative_v0.fits` (HDU QSO_CAT) lists every observation of every quasar target across surveys and
programs, so a TARGETID can appear several times; we keep ONE row per TARGETID, preferring main-survey dark-time
observations, then any main-survey row, then the rest. The sample is 1.9 < z < 3.1 so that pairs within
|r_par| < 30 Mpc/h cover the forest 1.96 <= z <= 3.0. Positions are in degrees; `qid` is the TARGETID, which
matches the forest LOS_ID of `desi_io.read_deltas` (used to drop the quasar's own sightline from the pairs).
"""
from __future__ import annotations
from dataclasses import dataclass
import sys
from pathlib import Path
import numpy as np
import fitsio
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline'):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import RAW
from cosmo import chi as chi_of_z

QSO_CAT=RAW/'desi/qso_iron/QSO_cat_iron_cumulative_v0.fits'


@dataclass
class QuasarSet:
    qid: np.ndarray      # TARGETID
    ra: np.ndarray       # degrees
    dec: np.ndarray
    z: np.ndarray
    chi: np.ndarray      # Mpc/h
    attrs: dict

    @property
    def n(self): return len(self.qid)


def read_quasars(zmin=1.9,zmax=3.1,path=QSO_CAT,verbose=True):
    t=fitsio.FITS(str(path))['QSO_CAT']
    d=t.read(columns=['TARGETID','TARGET_RA','TARGET_DEC','Z','SURVEY','PROGRAM'])
    survey=np.char.strip(d['SURVEY'].astype(str)); program=np.char.strip(d['PROGRAM'].astype(str))
    rank=np.where((survey=='main')&(program=='dark'),0,np.where(survey=='main',1,2))
    order=np.lexsort((rank,d['TARGETID']))            # by TARGETID, best rank first
    tid=d['TARGETID'][order]; first=np.r_[True,tid[1:]!=tid[:-1]]; keep=order[first]
    d=d[keep]; z=np.asarray(d['Z'],float); m=(z>zmin)&(z<zmax)&np.isfinite(z)
    q=QuasarSet(np.asarray(d['TARGETID'][m],np.int64),np.asarray(d['TARGET_RA'][m],float),np.asarray(d['TARGET_DEC'][m],float),z[m],
                np.asarray(chi_of_z(z[m]),np.float32),{'zmin':zmin,'zmax':zmax,'source':str(path),'rows_total':int(t.get_nrows()),
                                                       'unique_targets':int(len(keep)),'main_dark':int(((survey[keep]=='main')&(program[keep]=='dark')&m).sum())})
    if verbose: print(f'quasars: {t.get_nrows()} rows, {len(keep)} unique targets, {q.n} with {zmin} < z < {zmax} ({q.attrs["main_dark"]} main/dark)',flush=True)
    return q


if __name__=='__main__':
    import h5py
    q=read_quasars()
    with h5py.File('/data/LyaLenser/stageb/dr1_lowz_v7/sightlines.h5') as g: qid=g['sightlines/qid'][()]
    print('with a forest in the iteration-11 sample:',int(np.isin(q.qid,qid).sum()),'; forests whose quasar is in the sample:',int(np.isin(qid,q.qid).sum()))
