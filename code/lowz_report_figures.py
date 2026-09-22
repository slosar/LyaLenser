"""Figures of the dedicated low-redshift-tracer report (report/lowz/lowz.tex). Standalone; reads the committed JSON
products (report/lowz_validation.json, report/stageb/*.json), the tracer summary ($LYALENSER_DATA/lowz_split) and
the DR1 correlation table ($LYALENSER_DATA/stageb/dr1_lowz/xi.h5); writes report/lowz/figures/*.pdf."""
from pathlib import Path
import json, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
for p in (HERE,HERE/'pipeline'):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from cosmo import chi as chi_of_z
from cross_spectrum import kernel
OUT=ROOT/'report/lowz/figures'; OUT.mkdir(parents=True,exist_ok=True)
# Okabe-Ito, fixed order (colour-vision safe)
C={'blue':'#0072B2','orange':'#E69F00','green':'#009E73','red':'#D55E00','purple':'#CC79A7','sky':'#56B4E9','yellow':'#F0E442','black':'#000000'}
TRACER_COLOR={'LRG':C['red'],'ELG':C['blue'],'QSO':C['green'],'BGS':C['orange'],'BOSS':C['purple']}
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25,'grid.linewidth':.5,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})

LOWZ=next(q for q in (DATA/'lowz_v3',DATA/'lowz_v2',DATA/'lowz_split') if (q/'summary.json').exists())
lowz=json.load(open(LOWZ/'summary.json'))
cmb_planck=json.load(open(ROOT/'report/stageb/cmb_bias_check_planck.json')) if (ROOT/'report/stageb/cmb_bias_check_planck.json').exists() else None
cmb=json.load(open(ROOT/'report/stageb/cmb_bias_check.json'))
DR1JSON=next(q for q in (ROOT/'report/stageb/dr1_lowz_v7.json',ROOT/'report/stageb/dr1_lowz_v6.json',ROOT/'report/stageb/dr1_lowz_v4.json',ROOT/'report/stageb/dr1_lowz_v3.json',ROOT/'report/stageb/dr1_lowz.json') if q.exists())
dr1=json.load(open(DR1JSON))
mock=json.load(open(ROOT/'report/lowz_validation.json'))
tracers=[(sl,lab,t) for sl in lowz['slices'] for lab,t in sl['tracers'].items()]
order=sorted(range(len(tracers)),key=lambda i:(tracers[i][1].split('_')[0],tracers[i][0]['zmin']))
tracers=[tracers[i] for i in order]

# 1. tracer auto-spectra (cross of halves) with the fitted bias and the empirical shot noise
NROW=(len(tracers)+3)//4
fig,axes=plt.subplots(NROW,4,figsize=(11,2.5*NROW),sharex=True)
for ax,(sl,lab,t) in zip(axes.ravel(),tracers):
    L=np.array(t['L']); cx=np.array(t['cross_cl_binned']); T=np.array(t['theory_unit_bias_binned']); shot=t['shot_s_unit_bias']; f=t['bias_fit']
    band=(L>=f['band'][0])&(L<=f['band'][1]); col=TRACER_COLOR[lab.split('_')[0]]
    ax.plot(L[L>0],f['b']**2*T[L>0],color=col,lw=2,label=f"$b^2 C_\\ell$, $b={f['b']:.2f}\\pm{f['sigma_b']:.2f}$")
    ax.plot(L[L>0],cx[L>0],'o',ms=4,color=col,mfc='white',mew=1.2,label='cross of two random halves')
    ax.axhline(shot,color='0.4',ls='--',lw=1,label='shot noise (half difference)')
    ax.axvspan(f['band'][0],f['band'][1],color=col,alpha=.08)
    ax.set(yscale='log',xlim=(20,1000),title=lab.replace('_',' ',1).replace('_','–'),ylim=(max(1e-11,min(cx[L>40].min(),shot)/3),max(cx[L>0].max(),shot)*3))
    ax.legend(fontsize=6.5,loc='lower left')
for ax in axes[-1]: ax.set_xlabel(r'$\ell$')
for ax in axes[:,0]: ax.set_ylabel(r'$C_\ell$ (unit-bias $\kappa$ units)')
fig.suptitle('DESI DR1 and BOSS DR12 tracers: kernel-weighted maps, bias from the large-scale auto-spectrum (shaded band: $40\\leq\\ell\\leq0.2\\chi$)',fontsize=10)
fig.tight_layout(); fig.savefig(OUT/'tracer_spectra.pdf'); plt.close(fig)

# 2. ACT kappa x tracer cross-spectra
NROW=(len(tracers)+3)//4
fig,axes=plt.subplots(NROW,4,figsize=(11,2.5*NROW),sharex=True)
for ax,(sl,lab,t) in zip(axes.ravel(),tracers):
    c=cmb['tracers'][lab]; L=np.array(c['L']); cx=np.array(c['cross']); T=np.array(c['theory_unit_bias']); col=TRACER_COLOR[lab.split('_')[0]]
    var=None
    ax.plot(L,c['b_cmb_cross']*T,color=col,lw=2,label=f"$b_{{\\kappa}}={c['b_cmb_cross']:.2f}\\pm{c['sigma_b_cmb_cross']:.2f}$")
    ax.plot(L,c['b_auto']*T,color='0.3',ls=':',lw=1.2,label=f"$b_{{\\rm auto}}={c['b_auto']:.2f}$")
    ax.plot(L,cx,'o',ms=4,color=col,mfc='white',mew=1.2,label=r'$\langle m\,\kappa_{\rm CMB}\rangle$')
    ax.axhline(0,color='0.6',lw=.8); ax.set(title=lab.replace('_',' ',1).replace('_','–'),xlim=(20,700)); ax.legend(fontsize=6.5)
for ax in axes[-1]: ax.set_xlabel(r'$\ell$')
for ax in axes[:,0]: ax.set_ylabel(r'$C_\ell^{m\kappa}$')
fig.suptitle('Cross-spectra of the tracer maps with the ACT DR6 lensing map (joint sky fraction 0.05–0.08)',fontsize=10)
fig.tight_layout(); fig.savefig(OUT/'cmb_cross.pdf'); plt.close(fig)

# 3. bias comparison
fig,ax=plt.subplots(figsize=(6,3.2)); x=np.arange(len(tracers))
for i,(sl,lab,t) in enumerate(tracers):
    col=TRACER_COLOR[lab.split('_')[0]]; c=cmb['tracers'][lab]
    ax.errorbar(i-.12,t['bias_fit']['b'],t['bias_fit']['sigma_b'],fmt='o',color=col,ms=6,capsize=2,label='auto-spectrum (split halves)' if i==0 else None)
    ax.errorbar(i+.12,c['b_cmb_cross'],c['sigma_b_cmb_cross'],fmt='s',color=col,mfc='white',ms=6,capsize=2,label=r'ACT $\kappa$ cross' if i==0 else None)
    if cmb_planck is not None and lab in cmb_planck['tracers']:
        cp=cmb_planck['tracers'][lab]; ax.errorbar(i+.3,cp['b_cmb_cross'],cp['sigma_b_cmb_cross'],fmt='^',color=col,mfc='white',ms=5,capsize=2,label=r'Planck $\kappa$ cross' if i==0 else None)
ax.set_xticks(x); ax.set_xticklabels([lab.replace('_','\n',1).replace('_','–') for _,lab,_ in tracers],fontsize=7.5); ax.set_ylabel('linear bias $b$'); ax.set_ylim(0,3.2); ax.legend(loc='upper left',fontsize=8)
fig.tight_layout(); fig.savefig(OUT/'biases.pdf'); plt.close(fig)

# 4. mock validation: ensemble of the combined amplitude and the per-slice responses
S=mock['sparse']; a0=np.array([d['fits']['A0_R1']['combined']['A'] for d in S]); a1=np.array([d['fits']['A1_R1']['combined']['A'] for d in S])
sn=sorted(k for k in S[0]['fits']['A1_R1'] if k.startswith('slice'))
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9,3.4))
bins=np.linspace(-6,7,40)
ax1.hist(a0,bins,color=C['blue'],alpha=.55,label=f'$A_{{\\rm true}}=0$: ${a0.mean():.2f}\\pm{a0.std(ddof=1)/20:.2f}$')
ax1.hist(a1,bins,color=C['orange'],alpha=.55,label=f'$A_{{\\rm true}}=1$: ${a1.mean():.2f}\\pm{a1.std(ddof=1)/20:.2f}$')
ax1.axvline(0,color=C['blue'],lw=1); ax1.axvline(1,color=C['orange'],lw=1); ax1.set(xlabel='fitted $A$ (combined template)',ylabel='seeds (of 400)'); ax1.legend(fontsize=8)
labels=[f"{s['zmin']:g}–{s['zmax']:g}" for s in mock['sparse'][0]['templates']['slices']]
resp=[np.array([d['fits']['A1_R1'][n]['A']-d['fits']['A0_R1'][n]['A'] for d in S]) for n in sn]
ax2.errorbar(np.arange(len(sn)),[r.mean() for r in resp],[r.std(ddof=1)/20 for r in resp],fmt='o',color=C['green'],capsize=3,label='per slice')
ax2.errorbar([len(sn)],[(a1-a0).mean()],[(a1-a0).std(ddof=1)/20],fmt='D',color=C['red'],capsize=3,label='combined')
t1=np.array([d['fits']['A1_R1']['truth']['A']-d['fits']['A0_R1']['truth']['A'] for d in S]); ax2.axhline(t1.mean(),color='0.4',ls='--',lw=1,label=f'truth template {t1.mean():.3f}')
ax2.axhline(1,color='0.7',lw=.8); ax2.set_xticks(range(len(sn)+1)); ax2.set_xticklabels(labels+['all'],fontsize=8); ax2.set(ylabel='paired response $A(1)-A(0)$',ylim=(0.2,1.5)); ax2.legend(fontsize=8)
fig.suptitle('Mock validation (400 realisations, 400 deg$^2$, lognormal tracers of the same realisation)',fontsize=10)
fig.tight_layout(); fig.savefig(OUT/'mock_validation.pdf'); plt.close(fig)

# 5. DR1 result
fits=dr1['fits']; names=list(fits); rnd=np.array(dr1['random_templates']['amplitudes'])
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9,3.4),gridspec_kw={'width_ratios':[1.3,1]})
for i,n in enumerate(names):
    s=fits[n]; col=C['red'] if n=='combined' else C['blue']
    ax1.errorbar(i,s['A'],s['jk_error'],fmt='D' if n=='combined' else 'o',color=col,capsize=3,ms=6)
    ax1.errorbar(i,s['A'],s['sigma_F'],fmt='none',ecolor=col,elinewidth=3,alpha=.35)
ax1.errorbar(len(names),dr1['joint']['A'],dr1['joint']['error'],fmt='s',color=C['purple'],capsize=3,ms=6)
ax1.axhline(0,color='0.5',lw=.8); ax1.axhline(1,color='0.7',ls='--',lw=.8)
ax1.set_xticks(range(len(names)+1)); ax1.set_xticklabels(['combined']+[n.replace('slice_','').replace('_','–') for n in names[1:]]+['jackknife\ncombination'],fontsize=7.5,rotation=20)
ax1.set(ylabel='lensing amplitude $A$',ylim=(-5,9)); ax1.set_title('DR1, single slab $%g<z<%g$ (bars: jackknife; thick: Fisher)'%(dr1.get('zmin',2.1),dr1.get('zmax',3.0)),fontsize=9)
ax2.hist(rnd,np.linspace(-2.2,2.2,23),color='0.6',label=f"{len(rnd)} random templates: ${rnd.mean():.2f}\\pm{rnd.std(ddof=1)/np.sqrt(len(rnd)):.2f}$, scatter {rnd.std(ddof=1):.2f}")
ax2.axvline(fits['combined']['A'],color=C['red'],lw=2,label=f"data $A={fits['combined']['A']:.2f}\\pm{fits['combined']['jk_error']:.2f}$")
ax2.set(xlabel='$A$',ylabel='realisations'); ax2.legend(fontsize=7.5,loc='upper left')
fig.tight_layout(); fig.savefig(OUT/'dr1_result.pdf'); plt.close(fig)

# 6. DR1 forest correlation: the measured 1 Mpc/h table against the two-parameter fit and the corrected fit
from scipy.interpolate import RegularGridInterpolator
def read_table(f,group):
    g=f[group]; xi=g['xi'][()]; xirp=g['xi_rp'][()]
    if xi.ndim==3:      # iteration-10 layered table: show the layer at the reference distance chi(2.4)
        nodes=g['chi_nodes'][()]; t=np.clip((float(chi_of_z(2.4))-nodes[0])/(nodes[1]-nodes[0]),0,len(nodes)-1); i=min(int(t),len(nodes)-2); fr=t-i   # the layer at z = 2.4, the reference of the power laws
        xi=(1-fr)*xi[i]+fr*xi[i+1]; xirp=(1-fr)*xirp[i]+fr*xirp[i+1]
    return dict(rp=g['r_perp'][()],rz=g['r_par'][()],xi=xi,xirp=xirp,meta=json.loads(g.attrs['meta']))
STAGEB=DATA/'stageb'
RUN=next(q for q in (STAGEB/'dr1_lowz_v7',STAGEB/'dr1_lowz_v6',STAGEB/'dr1_lowz_v4',STAGEB/'dr1_lowz_v3',STAGEB/'dr1_lowz') if (q/'xi.h5').exists())
with h5py.File(RUN/'xi.h5') as f:
    T=read_table(f,'xi'); num=f['xi/coarse_num'][()]; den=f['xi/coarse_den'][()]
    T0=read_table(f,'xi_uncorrected') if 'xi_uncorrected' in f else None
TPREV=None
if RUN.name!='dr1_lowz' and (STAGEB/'dr1_lowz/xi.h5').exists():
    with h5py.File(STAGEB/'dr1_lowz/xi.h5') as f: TPREV=read_table(f,'xi')   # the iteration-5 table of the first run
raw=np.divide(num,den,out=np.zeros_like(num),where=den>0)
err=np.divide(1.,np.sqrt(den),out=np.full_like(den,np.inf),where=den>0)
centres=np.arange(raw.shape[0])+.5
CELL=tuple(np.meshgrid(np.arange(30)+.5,np.arange(30)+.5,indexing='ij'))
def cells(t,key='xi'):
    return RegularGridInterpolator((t['rp'],t['rz']),t[key],bounds_error=False,fill_value=np.nan)(CELL)
# With the redshift-evolving table (iteration 10+) the collapsed cells are a pair-weighted average over redshift,
# so the model to compare with is the same average of the layered table: the z-resolved cells are measured once
# from the run's sightlines (cached in xi.h5 'cells'), and the model layer at each cell's mean distance is
# coarse-binned and averaged with the cell weights.
def zcells(run):
    import h5py as _h5
    with _h5.File(run/'xi.h5','a') as f:
        if 'cells' in f and 'chisum' in f['cells']:
            g=f['cells']; return g['num'][()],g['den'][()],g['chisum'][()]
        from mock import load_sightlines
        from xi_model import xi_from_data
        from campaign4 import campaign_config
        fit=json.loads(f['xi'].attrs['meta'])['fit']; edges=tuple(fit['z_edges'])
        cfg=campaign_config(1.).copy(xi_z_evolution=True,xi_z_edges=edges,slabs=((edges[0],edges[-1]),))
        sl=load_sightlines(run/'sightlines.h5'); m=xi_from_data(sl,cfg); n_,d_,c_=m.counts_z['all']
        g=f.require_group('cells'); [g.__delitem__(k) for k in list(g)]; g['num']=n_; g['den']=d_; g['chisum']=c_; g['z_edges']=np.asarray(edges)
        return n_,d_,c_
def layered_cells(group,key='xi'):
    """Coarse cells [n_z, 30, 30] of a layered table at each z-bin cell's mean distance, and their den-weighted collapse."""
    with h5py.File(RUN/'xi.h5') as f:
        g=f[group]; arr=g[key][()]; nodes=g['chi_nodes'][()] if 'chi_nodes' in g else None; rp=g['r_perp'][()]; rz=g['r_par'][()]
    if nodes is None: return None
    NZ,DEN,CS=zcells(RUN); nz=NZ.shape[0]; out=np.zeros((nz,30,30))
    for k in range(nz):
        d=DEN[k][:30,:30]; cm=np.divide(CS[k][:30,:30],d,out=np.full_like(d,np.nan),where=d>0); cm=np.where(np.isfinite(cm),cm,np.nanmean(cm))
        t=np.clip((cm-nodes[0])/(nodes[1]-nodes[0]),0,len(nodes)-1); i=np.minimum(t.astype(int),len(nodes)-2); fr=t-i
        lay=lambda idx: RegularGridInterpolator((rp,rz),arr[idx],bounds_error=False,fill_value=np.nan)(CELL)
        # interpolate the layer per cell: evaluate every needed layer once
        vals=np.zeros((30,30))
        for j in np.unique(np.r_[i.ravel(),i.ravel()+1]):
            L=lay(j); w=np.where(i==j,1-fr,0)+np.where(i+1==j,fr,0); vals+=w*L
        out[k]=vals
    W=DEN[:,:30,:30]; return out,np.divide((out*W).sum(axis=0),W.sum(axis=0),out=np.full((30,30),np.nan),where=W.sum(axis=0)>0)
ZMODEL={}
for grp in ('xi','xi_uncorrected'):
    try:
        r_=layered_cells(grp)
        if r_ is not None: ZMODEL[grp]=r_[1]
    except Exception as e_: print('z-resolved model unavailable for',grp,e_)
def model_cells(t,grp):
    return ZMODEL[grp] if grp in ZMODEL else cells(t)
fig,axes=plt.subplots(1,3,figsize=(12,3.6),gridspec_kw={'width_ratios':[1.25,1,1]})
for (lo,hi),col in zip(((0,2),(4,6),(10,12)),(C['blue'],C['orange'],C['green'])):
    m=(centres>=lo)&(centres<hi); ok=(centres>=3)&(centres<30)
    y=raw[:,m].mean(axis=1); e=np.sqrt((err[:,m]**2).sum(axis=1))/m.sum()
    axes[0].errorbar(centres[ok],centres[ok]**2*y[ok],centres[ok]**2*e[ok],fmt='o',ms=3.5,color=col,capsize=0,
                     lw=1,label=f'$r_\\parallel$ {lo}--{hi}')
    for t,grp,ls,lw in ((T,'xi','-',1.6),(T0,'xi_uncorrected',':',1.2)):
        if t is None: continue
        if grp in ZMODEL:
            m30=m[:30]; ok30=ok[:30]; c30=centres[:30]; mc=ZMODEL[grp][:,m30].mean(axis=1); axes[0].plot(c30[ok30],c30[ok30]**2*mc[ok30],ls,color=col,lw=lw); continue
        mz=(t['rz']>=lo)&(t['rz']<hi); yf=t['xi'][:,mz].mean(axis=1); okf=(t['rp']>=3)&(t['rp']<30)
        axes[0].plot(t['rp'][okf],t['rp'][okf]**2*yf[okf],ls,color=col,lw=lw)
axes[0].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$r_\perp^2\,\xi_F$',xlim=(2,31))
axes[0].plot([],[],'-',color='0.3',lw=1.6,label='corrected table')
axes[0].plot([],[],':',color='0.3',lw=1.2,label='base model')
axes[0].legend(fontsize=7,ncol=2,loc='upper left'); axes[0].set_title('DR1 forest correlation',fontsize=9)
for ax,t,grp,ttl in ((axes[1],T0,'xi_uncorrected','base model (evolving Kaiser)'),(axes[2],T,'xi','with the spline correction')):
    if t is None: continue
    r=(raw[:30,:30]-model_cells(t,grp))/err[:30,:30]; r[:3]=np.nan
    im=ax.imshow(r.T,origin='lower',extent=(0,30,0,30),vmin=-5,vmax=5,cmap='RdBu_r')
    ax.set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$r_\parallel$ ($h^{-1}$Mpc)'); ax.grid(False)
    ax.set_title(f'residual / $\\sigma$, {ttl}',fontsize=9); fig.colorbar(im,ax=ax,shrink=.85)
fig.tight_layout(); fig.savefig(OUT/'dr1_xi.pdf'); plt.close(fig)

# 6b. what the correction does: the same-wavelength term and the change in the response kernel
REF=TPREV if TPREV is not None else T0
if T0 is not None and REF is not None:
    fit=T['meta']['fit']
    prof=np.array(json.load(open(ROOT/'report/signal_profile.json'))['information_2d']); prof=prof/prof.sum()
    gN=cells(T,'xirp'); gO=cells(REF,'xirp')
    W=den[:30,:30]*((np.arange(30)+.5)**2)[:,None]; W[:3]=0
    Afac=1/(np.sum(W*gO*gN)/np.sum(W*gO*gO))
    fig,axes=plt.subplots(1,3,figsize=(12,3.4))
    rp=np.arange(3,30)+.5
    f0=RegularGridInterpolator((T0['rp'],T0['rz']),T0['xi'])
    axes[0].errorbar(rp,1e3*(raw[3:30,0]-f0(np.column_stack((rp,np.full_like(rp,.5))))),1e3*err[3:30,0],
                     fmt='o',ms=3.5,color=C['blue'],capsize=0,lw=1,label='measured first radial bin\nminus the two-parameter fit')
    if 'same_wavelength' in fit:
        axes[0].plot(fit['same_wavelength_rperp'],1e3*np.array(fit['same_wavelength']),'-',color=C['red'],lw=1.6,
                     label='fitted same-wavelength term $N(r_\\perp)$')
    axes[0].axhline(0,color='0.6',lw=.8); axes[0].legend(fontsize=7)
    axes[0].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$10^3\,\Delta\xi_F$',
                title='same-wavelength excess ($r_\\parallel<1$)')
    ratio=np.where(np.abs(gO)>.15*np.max(np.abs(gO),axis=1,keepdims=True),gN/gO,np.nan); ratio[:3]=np.nan
    im=axes[1].imshow(ratio.T,origin='lower',extent=(0,30,0,30),vmin=.85,vmax=1.15,cmap='PuOr_r')
    axes[1].contour(np.arange(30)+.5,np.arange(30)+.5,prof.T,levels=[2e-3,8e-3],colors='k',linewidths=.7)
    axes[1].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r'$r_\parallel$ ($h^{-1}$Mpc)',ylim=(0,12),xlim=(3,30),
                title=r"kernel ratio $\xi'_{\rm new}/\xi'_{\rm published}$"); axes[1].grid(False)
    fig.colorbar(im,ax=axes[1],shrink=.85)
    gn=RegularGridInterpolator((T['rp'],T['rz']),T['xirp']); go=RegularGridInterpolator((REF['rp'],REF['rz']),REF['xirp'])
    x=np.linspace(3,30,80)
    for z,col in zip((0.5,3.5),(C['blue'],C['green'])):
        axes[2].plot(x,-1e3*x**2*go(np.column_stack((x,np.full_like(x,z)))),':',color=col,lw=1.3)
        axes[2].plot(x,-1e3*x**2*gn(np.column_stack((x,np.full_like(x,z)))),'-',color=col,lw=1.7,
                     label=f'$r_\\parallel={z:g}$')
    axes[2].plot([],[],'-',color='0.3',lw=1.7,label='corrected'); axes[2].plot([],[],':',color='0.3',lw=1.3,label='published')
    axes[2].set(xlabel=r'$r_\perp$ ($h^{-1}$Mpc)',ylabel=r"$-10^3 r_\perp^2\,\partial\xi_F/\partial r_\perp$",
                title=f'response kernel (amplitude $\\times${Afac:.3f})'); axes[2].legend(fontsize=7,ncol=2)
    fig.tight_layout(); fig.savefig(OUT/'xi_correction.pdf'); plt.close(fig)
    print('xi fit: chi2 %.0f -> %.0f over %d cells; kernel A factor %.4f'%(T0['meta']['fit']['chi2'],fit['chi2'],fit['cells'],Afac))

# 7. lensing kernel and the slices
from lowz import slice_spectra
# the slices of the tracer templates actually used (six with BGS and BOSS; five in the first version)
sl=[{'zmin':s['zmin'],'zmax':s['zmax'],'tracers':list(s['tracers'])} for s in lowz['slices']]
ZREF=float(lowz.get('z_ref',2.4)); cref=float(chi_of_z(ZREF)); cc=np.linspace(50,cref-1,600); W=kernel(cc,cref)
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9,3.2))
ax1.plot(cc,W/W.max(),color='k',lw=1.5,label=r'$W_{\rm Ly\alpha}(\chi)$, source at $z=%.2f$'%ZREF)
zt=[0.1,0.4,0.8,1.1,1.6,2.1]; ax1.set_xticks([float(chi_of_z(z)) for z in zt]); ax1.set_xticklabels([f'z={z:g}' for z in zt],fontsize=7.5)
cols=[C['yellow'],C['blue'],C['orange'],C['green'],C['purple'],C['sky']][-len(sl):]; frac=[]; Lq=[40,100,300]
for s,col in zip(sl,cols):
    c1,c2=float(chi_of_z(s['zmin'])),float(chi_of_z(s['zmax'])); ax1.axvspan(c1,c2,color=col,alpha=.3)
    ax1.text(.5*(c1+c2),1.03,'+'.join(sorted({t.split('_')[0] for t in s['tracers']})),ha='center',va='bottom',fontsize=6.5,rotation=0)
    L,Cs=slice_spectra(s['zmin'],s['zmax'],cref,400); ll=np.interp(Lq,L,Cs[:,1,1]); frac.append(ll)
frac=np.array(frac)
from three_tracer import limber
wl=lambda c: kernel(c,cref); full=limber(np.array(Lq,float),wl,wl,1.,cref-1,nchi=2000,to_recombination=False)
ax1.set(xlabel=r'comoving distance $\chi$ ($h^{-1}$Mpc)',ylabel='normalised kernel',ylim=(0,1.2)); ax1.legend(fontsize=8,loc='upper left')
n=len(sl)
for k,(Lv,col) in enumerate(zip(Lq,(C['blue'],C['orange'],C['green']))):
    ax2.bar(np.arange(n)+(k-1)*.25,frac[:,k]/full[k],width=.25,color=col,label=f'$L={Lv}$: total {frac[:,k].sum()/full[k]:.2f}')
ax2.set_xticks(range(n)); ax2.set_xticklabels([f"{s['zmin']:g}–{s['zmax']:g}" for s in sl],fontsize=8); ax2.set(ylabel=r'fraction of $C_L^{\kappa\kappa}$ of $\kappa_{\rm Ly\alpha}$',xlabel='tracer slice'); ax2.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT/'kernel_slices.pdf'); plt.close(fig)
print('fractions of kappa_lya power per slice at L=40,100,300:',np.round(frac/full,3).tolist(),'total',np.round((frac/full).sum(axis=0),3))
import shutil
for f in ('signal_profile.pdf','lowz_iteration7_normalisation.pdf'): shutil.copy(ROOT/'report/figures'/f,OUT/f)
print('figures written to',OUT)
