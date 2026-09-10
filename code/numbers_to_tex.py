"""Convert report/numbers.json (written by make_plots.py) into LaTeX macros and a table body."""
import json, numpy as np
d = json.load(open("../report/numbers.json"))
def sci(x, digits=2):
    if x is None or not np.isfinite(x): return r"\infty"
    m, e = f"{x:.{digits}e}".split("e"); return rf"{m}\times10^{{{int(e)}}}"
def fl(x, digits=1):
    return "--" if x is None else f"{x:.{digits}f}"
with open("../report/numbers.tex", "w") as f:
    f.write(f"\\newcommand{{\\chiStar}}{{{d['chi']:.0f}}}\n")
    f.write(f"\\newcommand{{\\NcAct}}{{{sci(d['Nc_act'])}}}\n\\newcommand{{\\NcPl}}{{{sci(d['Nc_pl'])}}}\n")
    f.write(f"\\newcommand{{\\Cll}}{{{sci(d['Cll100'])}}}\n\\newcommand{{\\Clc}}{{{sci(d['Clc100'])}}}\n\\newcommand{{\\Ccc}}{{{sci(d['Ccc100'])}}}\n")
    f.write(f"\\newcommand{{\\lmaxFifty}}{{{d['lmax']['50']:.0f}}}\n")
    f.write(f"\\newcommand{{\\NFifty}}{{{sci(d['N100']['50_2.0'])}}}\n\\newcommand{{\\NbhFifty}}{{{sci(d['Nbh100']['50_2.0'])}}}\n")
    f.write(f"\\newcommand{{\\snrFiftyAct}}{{{fl(d['snr']['50_N_ACT'],0)}}}\n\\newcommand{{\\snrFiftyPl}}{{{fl(d['snr']['50_N_Planck'],0)}}}\n")
    f.write(f"\\newcommand{{\\snrFiftyBhAct}}{{{fl(d['snr']['50_N_bh_slice_ACT'],0)}}}\n")
with open("../report/numbers_table.tex", "w") as f:
    f.write("\\begin{tabular}{rrrccrrr}\n\\toprule\n$n_q$ [deg$^{-2}$] & $l_{\\max}$ & $k_{\\parallel,\\max}$ & $N_\\kappa(100)$ & $N_\\kappa^{\\rm BH}(100)$ & S/N ACT & S/N Planck & S/N ACT, BH\\\\\n\\midrule\n")
    for n in ["20", "50", "100", "400"]:
        for km in ["1.0", "2.0", "5.0"]:
            key = f"{n}_{km}"
            row = [n if km == "1.0" else "", f"{d['lmax'][n]:.0f}" if km == "1.0" else "", km,
                   "$" + sci(d["N100"][key]) + "$", "$" + sci(d["Nbh100"][key]) + "$"]
            if km == "2.0":
                row += [fl(d['snr'][f'{n}_N_ACT']), fl(d['snr'][f'{n}_N_Planck']), fl(d['snr'][f'{n}_N_bh_slice_ACT'])]
            else:
                row += ["", "", ""]
            f.write(" & ".join(row) + r" \\" + "\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
print(open("../report/numbers.tex").read())
