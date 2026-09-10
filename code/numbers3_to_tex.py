"""Convert report/numbers3.json (from three_tracer.py) into LaTeX macros and a table."""
import json, numpy as np
d = json.load(open("../report/numbers3.json"))
def sci(x, digits=2):
    if x is None or not np.isfinite(x): return r"\infty"
    m, e = f"{x:.{digits}e}".split("e"); return rf"{m}\times10^{{{int(e)}}}"
s100 = d["spec"]["100"]; s40 = d["spec"]["40"]; p = d["params"]
with open("../report/numbers3.tex", "w") as f:
    f.write(f"\\newcommand{{\\Dslab}}{{{p['D_SLAB']:.0f}}}\n")
    f.write(f"\\newcommand{{\\betaHundred}}{{{s100['beta']:.3f}}}\n\\newcommand{{\\betaForty}}{{{s40['beta']:.3f}}}\n")
    f.write(f"\\newcommand{{\\rtwoQK}}{{{s100['r2_qk']:.3f}}}\n\\newcommand{{\\rtwoQKforty}}{{{s40['r2_qk']:.3f}}}\n")
    f.write(f"\\newcommand{{\\respOverLensForty}}{{{s40['resp_over_lens']:.2f}}}\n\\newcommand{{\\respOverLensHundred}}{{{s100['resp_over_lens']:.2f}}}\n")
    f.write(f"\\newcommand{{\\dbetaMag}}{{{100*s100['dbeta_over_beta']:.0f}}}\n")
    f.write(f"\\newcommand{{\\betaShot}}{{{sci(s100['beta2_shot'])}}}\n\\newcommand{{\\CkckcHundred}}{{{sci(s100['kckc'])}}}\n")
    f.write(f"\\newcommand{{\\CdkcHundred}}{{{sci(s100['dkc'])}}}\n\\newcommand{{\\CddHundred}}{{{sci(s100['dd'])}}}\n")
    rows = {(r["neff"], r["pn"], r.get("sigma_ln", 0.0)): r for r in d["table"]}
    for sig, tag in ((1.0, "One"), (2.0, "Two")):
        rw = rows[(25.0, 0.33, sig)]; rw50 = rows[(50.0, 0.33, sig)]
        f.write(f"\\newcommand{{\\snrDepW{tag}Act}}{{{rw['snr_deproj_ACT']:.1f}}}\n\\newcommand{{\\NW{tag}}}{{{sci(rw['N100'])}}}\n")
        f.write(f"\\newcommand{{\\snrDepWFifty{tag}Act}}{{{rw50['snr_deproj_ACT']:.1f}}}\n")
        nf = d["neff_factor"][str(sig)]
        f.write(f"\\newcommand{{\\nfac{tag}Low}}{{{nf['0.1']:.1f}}}\n\\newcommand{{\\nfac{tag}High}}{{{nf['1.0']:.1f}}}\n")
    r0 = rows[(25.0, 0.0, 0.0)]; r1 = rows[(25.0, 0.33, 0.0)]; r2 = rows[(25.0, 0.17, 0.0)]
    f.write(f"\\newcommand{{\\NnoiselessTF}}{{{sci(r0['N100'])}}}\n\\newcommand{{\\NnoisyTF}}{{{sci(r1['N100'])}}}\n")
    f.write(f"\\newcommand{{\\biasRatioTF}}{{{r1['bias_over_signal100']:.1f}}}\n")
    for tag, r in (("Noisy", r1), ("Clean", r2), ("Noiseless", r0)):
        f.write(f"\\newcommand{{\\snrDep{tag}Act}}{{{r['snr_deproj_ACT']:.1f}}}\n\\newcommand{{\\snrDep{tag}Pl}}{{{r['snr_deproj_Planck']:.1f}}}\n")
        f.write(f"\\newcommand{{\\snrBh{tag}Act}}{{{r['snr_bh_slice_ACT']:.1f}}}\n\\newcommand{{\\snrNaive{tag}Act}}{{{r['snr_naive_ACT']:.1f}}}\n")
    rf = rows[(50.0, 0.17, 0.0)]
    f.write(f"\\newcommand{{\\snrDepFutureAct}}{{{rf['snr_deproj_ACT']:.1f}}}\n")
with open("../report/numbers3_table.tex", "w") as f:
    f.write("\\footnotesize\\setlength{\\tabcolsep}{3.5pt}\n\\begin{tabular}{rrrcccrrrrr}\n\\toprule\n")
    f.write("$n_{\\rm eff}$ & $P_N$ & $\\sigma_{\\ln}$ & $N_\\kappa(100)$ & $N_\\kappa^{\\rm BH}(100)$ & bias/signal & \\multicolumn{5}{c}{cumulative S/N (ACT unless noted)}\\\\\n")
    f.write("{}[deg$^{-2}$] & [$h^{-1}$Mpc] & & & & at $L=100$ & naive & BH sl. & BH gl. & deproj. & dep. Pl.\\\\\n\\midrule\n")
    for r in d["table"]:
        f.write(f"{r['neff']:.0f} & {r['pn']:.2f} & {r.get('sigma_ln', 0.0):.0f} & ${sci(r['N100'])}$ & ${sci(r['Nbh100'])}$ & {r['bias_over_signal100']:.2f} & "
                f"{r['snr_naive_ACT']:.1f} & {r['snr_bh_slice_ACT']:.1f} & {r['snr_bh_global_ACT']:.1f} & {r['snr_deproj_ACT']:.1f} & {r['snr_deproj_Planck']:.1f}\\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print("ok")
