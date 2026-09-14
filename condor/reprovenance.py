"""Rewrite the stored provenance of campaign products after a change of the fingerprint FUNCTION only.

Used once (2026-09-12, iteration 4) when the machine-specific path fields were removed from the fingerprint so
that products from RACF and Perlmutter could be merged. It re-stamps every complete.json / attempt.json under the
mock root with the current fingerprint, in dependency order (dev -> freeze -> seeds/controls), updating the chained
digests (freeze.development_markers, seed/control provenance['freeze']). It never touches any product data.

Review 5 (finding 7) asked for a narrow tool: it now (1) refuses to run without an explicit migration manifest
naming every allowed difference with its old and new value, (2) verifies every artifact hash of every marker
before anything is rewritten, and (3) keeps the original provenance under ``provenance_before_migration`` in each
marker together with the manifest path. Float config fields are only allowed to differ by the 12-significant-digit
rounding introduced in the fingerprint.

Manifest (JSON): {"sources": {"code/pipeline/campaign4.py": ["<old sha256>", "<new sha256>"]},
                  "config_removed": ["data_root", "report_root"], "note": "..."}
Usage: python condor/reprovenance.py --mock-root $LYALENSER_DATA/mocks/iteration4 --manifest m.json [--scale 1] [--dry-run]
"""
import argparse, json, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'code/pipeline')); sys.path.insert(0,str(HERE.parent/'code'))
import campaign4 as c

def differences(old,new):
    diff=[]
    for k in set(old)|set(new):
        if k in ('config','sources'):
            for kk in set(old.get(k,{}))|set(new.get(k,{})):
                if old.get(k,{}).get(kk)!=new.get(k,{}).get(kk): diff.append(f'{k}.{kk}')
        elif old.get(k)!=new.get(k): diff.append(k)
    return diff

def close(old,new,key):
    """Float config entries that only differ by the 12-significant-digit rounding introduced in the fingerprint."""
    a=old.get('config',{}).get(key); b=new.get('config',{}).get(key)
    return isinstance(a,float) and isinstance(b,float) and abs(a-b)<=1e-9*max(abs(a),abs(b),1e-300)

def check(old,new,label,manifest):
    """Every difference must be named in the manifest with exactly the stored (old) and current (new) values."""
    bad=[]
    for d in differences(old,new):
        if d.startswith('config.'):
            key=d[7:]
            if key in manifest.get('config_removed',[]) and key not in new.get('config',{}): continue
            if close(old,new,key): continue
            bad.append(d)
        elif d.startswith('sources.'):
            key=d[8:]; pair=manifest.get('sources',{}).get(key)
            if pair and old.get('sources',{}).get(key)==pair[0] and new.get('sources',{}).get(key)==pair[1]: continue
            bad.append(d)
        else: bad.append(d)
    if bad: raise SystemExit(f'{label}: fingerprint differs in keys not covered by the manifest {bad}; refusing')

def verify_artifacts(path):
    d=json.loads(path.read_text())
    for name,h in d.get('artifacts',{}).items():
        p=path.parent/name
        if not p.exists() or c.digest(p)!=h: raise SystemExit(f'{path}: artifact {name} missing or changed; refusing')
    return d

def restamp(path,new_prov,dry,manifest_path):
    d=json.loads(path.read_text())
    if dry: return d
    d.setdefault('provenance_before_migration',d['provenance']); d['migration_manifest']=str(manifest_path)
    d['provenance']=new_prov; c.v.dump(path,d); return d

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--mock-root',type=Path,required=True); ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--manifest',type=Path,required=True); ap.add_argument('--dry-run',action='store_true'); a=ap.parse_args()
    manifest=json.loads(a.manifest.read_text())
    root=a.mock_root; base=c.fingerprint(a.scale)
    markers=[p for p in [root/f'dev/{s}/complete.json' for s in c.DEV_SEEDS]+[root/'freeze/complete.json'] if p.exists()]
    markers+=[p for sub in ('sparse','dense','controls') for p in sorted(root.glob(f'{sub}/*/complete.json'))]
    # 0. verify everything before rewriting anything
    for p in markers: verify_artifacts(p)
    for p in markers:
        old=json.loads(p.read_text())['provenance']
        check({k:v for k,v in old.items() if k not in ('freeze','basis','variant','seed','control')},base,str(p),manifest)
    # 1. development seeds
    for s in c.DEV_SEEDS:
        p=root/f'dev/{s}/complete.json'
        if p.exists(): restamp(p,base,a.dry_run,a.manifest); print('restamped',p)
    # 2. freeze: provenance + development markers
    fz=root/'freeze/complete.json'
    if fz.exists():
        d=json.loads(fz.read_text())
        d.setdefault('provenance_before_migration',d['provenance']); d['migration_manifest']=str(a.manifest); d['provenance']=base
        d['result']['development_markers']={str(s):c.digest(root/f'dev/{s}/complete.json') for s in c.DEV_SEEDS if (root/f'dev/{s}/complete.json').exists()}
        if (root/'basis/complete.json').exists(): d['result']['basis_marker']=c.digest(root/'basis/complete.json')
        for k in ('development.json','frozen.json'):
            q=root/'freeze'/k
            if q.exists() and not a.dry_run:
                if k=='frozen.json': c.v.dump(q,base)
                else:
                    dj=json.loads(q.read_text()); dj['development_markers']=d['result']['development_markers']; c.v.dump(q,dj)
        if not a.dry_run:
            d['artifacts']={str(p.relative_to(root/'freeze')):c.digest(p) for p in sorted((root/'freeze').rglob('*')) if p.is_file() and p.name!='complete.json'}
            c.v.dump(fz,d)
        print('restamped',fz)
    freeze_digest=c.digest(fz) if fz.exists() else None
    basis_digest=c.digest(root/'basis/complete.json') if (root/'basis/complete.json').exists() else None
    # 3. seeds and controls
    for sub in ('sparse','dense','controls'):
        for p in sorted(root.glob(f'{sub}/*/complete.json')):
            old=json.loads(p.read_text())['provenance']
            extra={k:old[k] for k in ('variant','seed','control') if k in old}
            new=base|{'freeze':freeze_digest}|({'basis':basis_digest} if basis_digest else {})|extra
            restamp(p,new,a.dry_run,a.manifest); print('restamped',p)
        for p in sorted(root.glob(f'{sub}/*/attempt.json')):
            if not (p.parent/'complete.json').exists():
                if not a.dry_run: p.unlink()
                print('removed stale attempt',p)
    print('DRY RUN, nothing written' if a.dry_run else 'done')

if __name__=='__main__': main()
