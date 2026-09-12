"""Rewrite the stored provenance of iteration-4 campaign products after a change of the fingerprint FUNCTION.

Used once (2026-09-12) when the machine-specific path fields were removed from the fingerprint so that products
from RACF and Perlmutter can be merged. It re-stamps every complete.json / attempt.json under the mock root with
the current fingerprint, in dependency order (dev -> freeze -> seeds/controls), updating the chained digests
(freeze.development_markers, seed/control provenance['freeze']). It never touches any product data, and it
refuses to run unless every difference between the stored and the current fingerprint is confined to the keys
listed in ALLOWED (the config path fields and the hashes of the files named in ALLOWED_SOURCES).

Usage: python condor/reprovenance.py --mock-root $LYALENSER_DATA/mocks/iteration4 [--scale 1] [--dry-run]
"""
import argparse, json, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'code/pipeline')); sys.path.insert(0,str(HERE.parent/'code'))
import campaign4 as c

ALLOWED_CONFIG={'data_root','report_root'}
ALLOWED_SOURCES={'code/pipeline/campaign4.py'}

def differences(old,new):
    diff=[]
    for k in set(old)|set(new):
        if k=='config':
            for kk in set(old.get('config',{}))|set(new.get('config',{})):
                if old.get('config',{}).get(kk)!=new.get('config',{}).get(kk): diff.append(f'config.{kk}')
        elif k=='sources':
            for kk in set(old.get('sources',{}))|set(new.get('sources',{})):
                if old.get('sources',{}).get(kk)!=new.get('sources',{}).get(kk): diff.append(f'sources.{kk}')
        elif old.get(k)!=new.get(k): diff.append(k)
    return diff

def close(old,new,key):
    """Float config entries that only differ by the 12-significant-digit rounding introduced in the fingerprint."""
    a=old.get('config',{}).get(key); b=new.get('config',{}).get(key)
    return isinstance(a,float) and isinstance(b,float) and abs(a-b)<=1e-9*max(abs(a),abs(b),1e-300)

def check(old,new,label):
    bad=[d for d in differences(old,new)
         if not (d.startswith('config.') and (d[7:] in ALLOWED_CONFIG or close(old,new,d[7:])))
         and not (d.startswith('sources.') and d[8:] in ALLOWED_SOURCES)]
    if bad: raise SystemExit(f'{label}: fingerprint differs in disallowed keys {bad}; refusing')

def restamp(path,new_prov,dry):
    d=json.loads(path.read_text())
    if dry: return d
    d['provenance']=new_prov; c.v.dump(path,d); return d

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--mock-root',type=Path,required=True); ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--dry-run',action='store_true'); a=ap.parse_args()
    root=a.mock_root; base=c.fingerprint(a.scale)
    # 1. development seeds
    for s in c.DEV_SEEDS:
        p=root/f'dev/{s}/complete.json'
        if p.exists():
            old=json.loads(p.read_text())['provenance']; check(old,base,str(p)); restamp(p,base,a.dry_run); print('restamped',p)
    # 2. freeze: provenance + development markers
    fz=root/'freeze/complete.json'
    if fz.exists():
        d=json.loads(fz.read_text()); check(d['provenance'],base,str(fz))
        d['provenance']=base
        d['result']['development_markers']={str(s):c.digest(root/f'dev/{s}/complete.json') for s in c.DEV_SEEDS if (root/f'dev/{s}/complete.json').exists()}
        for k in ('development.json','frozen.json'):
            q=root/'freeze'/k
            if q.exists():
                if k=='frozen.json' and not a.dry_run: c.v.dump(q,base)
                if k=='development.json' and not a.dry_run:
                    dj=json.loads(q.read_text()); dj['development_markers']=d['result']['development_markers']; c.v.dump(q,dj)
        # artifacts of the freeze directory changed (development.json, frozen.json) -> recompute their hashes
        if not a.dry_run:
            d['artifacts']={str(p.relative_to(root/'freeze')):c.digest(p) for p in sorted((root/'freeze').rglob('*')) if p.is_file() and p.name!='complete.json'}
            c.v.dump(fz,d)
        print('restamped',fz)
    freeze_digest=c.digest(fz) if fz.exists() else None
    # 3. seeds and controls
    for sub in ('sparse','dense','controls'):
        for p in sorted(root.glob(f'{sub}/*/complete.json')):
            d=json.loads(p.read_text()); old=d['provenance']
            extra={k:old[k] for k in ('variant','seed','control') if k in old}
            new=base|{'freeze':freeze_digest}|extra
            check({k:v for k,v in old.items() if k not in ('freeze',)},{k:v for k,v in new.items() if k not in ('freeze',)},str(p))
            restamp(p,new,a.dry_run); print('restamped',p)
        for p in sorted(root.glob(f'{sub}/*/attempt.json')):
            if not (p.parent/'complete.json').exists():
                # an interrupted phase: drop the stale attempt so that it restarts cleanly
                if not a.dry_run: p.unlink()
                print('removed stale attempt',p)
    print('DRY RUN, nothing written' if a.dry_run else 'done')

if __name__=='__main__': main()
