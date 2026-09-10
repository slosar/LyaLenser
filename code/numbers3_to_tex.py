"""Convert report/numbers3.json (from three_tracer.py) into LaTeX macros and a table (post review 1)."""
import json, numpy as np
d = json.load(open("../report/numbers3.json"))
def sci(x, digits=2):
    if x is None or not np.isfinite(x): return r"\infty"
    m, e = f"{x:.{digits}e}".split("e"); return rf"{m}\times10^{{{int(e)}}}"
s100, s40 = d["spec"]["100"], d["spec"]["40"]; p = d["params"]; mt = d["matched"]
rows = {(r["neff"], r["pn"], r.get("sigma_ln", 0.0)): r for r in d["table"]}
n0 = min(r["neff"] for r in d["table"])
r0, r1, r2 = rows[(n0, 0.0, 0.0)], rows[(n0, 0.33, 0.0)], rows[(n0, 0.17, 0.0)]
rf = rows[(50.0, 0.17, 0.0)]; rw = rows[(n0, 0.33, 2.0)]
with open("../report/numbers3.tex", "w") as f:
    W = lambda name, val: f.write(f"\\newcommand{{\\{name}}}{{{val}}}\n")
    W("Dslab", f"{p['D_SLAB']:.0f}"); W("nlos", f"{n0:.0f}")
    W("betaHundred", f"{s100['beta']:.3f}"); W("rtwoQK", f"{s100['r2_qk']:.3f}")
    W("respOverLensForty", f"{s40['resp_over_lens']:.2f}"); W("respOverLensHundred", f"{s100['resp_over_lens']:.2f}")
    W("dbetaMag", f"{100*s100['dbeta_over_beta']:.0f}"); W("betaShot", sci(s100['beta2_shot'])); W("shotS", sci(mt['shot_s']))
    W("CkckcHundred", sci(s100['kckc'])); W("CdkcHundred", sci(s100['dkc'])); W("CddHundred", sci(s100['dd']))
    W("Wlo", sci(mt['W_range'][0])); W("Whi", sci(mt['W_range'][1])); W("sigSub", f"{100*mt['skl_over_klkc100']:.1f}")
    W("NnoiselessTF", sci(r0['N100'])); W("NnoisyTF", sci(r1['N100'])); W("biasRatioTF", f"{r1['bias_over_signal100']:.2f}")
    for tag, r in (("Noisy", r1), ("Clean", r2), ("Noiseless", r0), ("Future", rf), ("WTwo", rw)):
        W(f"snrMat{tag}Act", f"{r['snr_matched_ACT']:.1f}"); W(f"snrMat{tag}Pl", f"{r['snr_matched_Planck']:.1f}")
        W(f"snrMatL{tag}Act", f"{r['snrmatched_ACT_L40_300']:.1f}"); W(f"snrMatL{tag}Pl", f"{r['snrmatched_Planck_L40_300']:.1f}")
        W(f"snrNaive{tag}Act", f"{r['snr_naive_ACT']:.1f}"); W(f"snrNaiveL{tag}Act", f"{r['snrnaive_ACT_L40_300']:.1f}")
        W(f"snrBh{tag}Act", f"{r['snr_bh_slice_ACT']:.1f}"); W(f"snrDep{tag}Act", f"{r['snr_deproj_ACT']:.1f}")
        W(f"snrRespL{tag}Act", f"{r['snrresp_ACT_L40_300']:.1f}"); W(f"snrResp{tag}Act", f"{r['snr_resp_ACT']:.1f}")
    for sig, tag in ((1.0, "One"), (2.0, "Two")):
        nf = d["neff_factor"][str(sig)]
        W(f"nfac{tag}Low", f"{nf['0.1']:.1f}"); W(f"nfac{tag}High", f"{nf['1.0']:.1f}")
with open("../report/numbers3_table.tex", "w") as f:
    f.write("\\footnotesize\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{rrrccrrrrrr}\n\\toprule\n")
    f.write("$\\bar n_{\\rm los}$ & $P_N$ & $\\sigma_{\\ln}$ & $N_\\kappa(100)$ & bias/signal & \\multicolumn{6}{c}{cumulative S/N}\\\\\n")
    f.write(" & & & & at $L=100$ & \\multicolumn{3}{c}{all $L$, ACT} & \\multicolumn{3}{c}{$40\\le L\\le300$}\\\\\n")
    f.write("{}[deg$^{-2}$] & [$h^{-1}$Mpc] & & & & naive & BH & matched & matched ACT & matched Pl. & response\\\\\n\\midrule\n")
    for r in d["table"]:
        f.write(f"{r['neff']:.0f} & {r['pn']:.2f} & {r.get('sigma_ln', 0.0):.0f} & ${sci(r['N100'])}$ & {r['bias_over_signal100']:.2f} & "
                f"{r['snr_naive_ACT']:.1f} & {r['snr_bh_slice_ACT']:.1f} & {r['snr_matched_ACT']:.1f} & "
                f"{r['snrmatched_ACT_L40_300']:.1f} & {r['snrmatched_Planck_L40_300']:.1f} & {r['snrresp_ACT_L40_300']:.1f}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print("ok")
