"""Result figures of the paper (Figs. 3-11) in Origin style (closed frame, inward ticks on all sides, Times New Roman),
sized to their printed width so that lettering prints at 8-10 pt."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullLocator

ROOT = Path(__file__).resolve().parent
RES, FIG = ROOT / "Results", ROOT / "Figures"
FIG.mkdir(exist_ok=True)
CM = 1 / 2.54
W2 = 14.5 * CM      # full text width
W1 = 9.0 * CM       # single-panel width

# Origin-like strong colours, one per wiring; variants differ by line style and symbol
C = {"2T": "#1f3fd8", "VMM": "#e0201b", "3T": "#008c3a", "4T": "#000000"}
STYLE = {  # key -> (label, colour, linestyle, marker)
    "2T": ("2T series", C["2T"], "-", "o"),
    "VMM": (r"VMM, $n^{*}$", C["VMM"], "-", "s"),
    "VMM-n2": (r"VMM, $n$ = 2", C["VMM"], (0, (5, 2)), "s"),
    "3T-L1": ("3T-L1", C["3T"], "-", "^"),
    "3T-L2": ("3T-L2a", C["3T"], (0, (5, 2)), "v"),
    "3T-L2b": ("3T-L2b", C["3T"], (0, (8, 2, 2, 2)), "<"),
    "3T-L3": ("3T-L3", C["3T"], (0, (1.5, 1.5)), ">"),
    "3T-L4": ("3T-L4", C["3T"], (0, (3, 1.5, 1, 1.5)), "D"),
    "4T": ("4T", C["4T"], "-", "o"),
}
SEQ_CMAP = LinearSegmentedColormap.from_list(
    "originblue", ["#ffffff", "#d7e6fa", "#9ec5f4", "#5598e7", "#256abf", "#163f7a", "#0b1f3d"])

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman"], "mathtext.fontset": "stix",
    "font.size": 9, "axes.labelsize": 9.5, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 7.5,
    "axes.linewidth": 1.0, "axes.edgecolor": "black", "axes.labelcolor": "black", "text.color": "black",
    "xtick.color": "black", "ytick.color": "black", "axes.grid": False,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 4, "ytick.major.size": 4, "xtick.minor.size": 2.2, "ytick.minor.size": 2.2,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9, "xtick.minor.width": 0.7, "ytick.minor.width": 0.7,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.spines.top": True, "axes.spines.right": True,
    "lines.linewidth": 1.3, "lines.markersize": 4, "lines.markeredgewidth": 0.6,
    "legend.frameon": False, "legend.handlelength": 2.6, "legend.borderaxespad": 0.4, "legend.labelspacing": 0.3,
    "savefig.dpi": 600, "figure.dpi": 150,
})


def read(name):
    with open(RES / f"{name}.csv", newline="") as f:
        rows = list(csv.reader(f))

    def conv(x):
        try:
            return float(x)
        except ValueError:
            return x
    return [dict(zip(rows[0], map(conv, r))) for r in rows[1:]]


def head(ax, s):
    """Panel label outside the frame: bold, above the top-left corner of the axes."""
    ax.set_title(s, loc="left", fontsize=10, fontweight="bold", pad=5)


def tag(ax, s, x=0.03, y=0.95, ha="left", va="top"):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", ha=ha, va=va)


def note(ax, s, x, y, ha="left", va="top", size=8):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=size, ha=ha, va=va)


def line(ax, x, y, key, **kw):
    lab, col, ls, mk = STYLE[key]
    opts = dict(color=col, ls=ls, marker=mk, mfc=col, mec=col, label=lab)
    opts.update(kw)
    return ax.plot(x, y, **opts)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", facecolor="white")
    plt.close(fig)


ref = json.load(open(RES / "reference_operating_points.json"))
NOPT, PIN = ref["nOpt"], ref["IncidentPower_W_m2"]
I = ref["Ideal"]
P0 = 100 * I["4T0"]["P"] / PIN

# ---------------------------------------------------------------- Fig: reference J-V
pe = read("reference_subcell_JV"); si = read("reference_si_JV")
t2 = read("reference_2T_JV"); vm = read("reference_VMM_JV")
fig, axs = plt.subplots(1, 2, figsize=(W2, 6.0 * CM), constrained_layout=True)
ax = axs[0]
ax.plot([r["V"] for r in pe], [r["J_top_mAcm2"] for r in pe], color=C["2T"], lw=1.6, label="QIBC perovskite top cell")
ax.plot([r["V"] for r in si], [r["J_bot_mAcm2"] for r in si], color=C["VMM"], lw=1.6, label="HIBC Si bottom cell (filtered)")
ax.plot([I["4T0"]["VTop"], I["4T0"]["VBot"]], [0.1 * I["4T0"]["JTop"], 0.1 * I["4T0"]["JBot"]], "o", color="black",
        mfc="white", ms=5.5, mew=1.1, label="Maximum power point")
ax.set(xlim=(0, 1.4), ylim=(0, 31), xlabel="Voltage (V)", ylabel=r"Current density (mA cm$^{-2}$)")
ax.legend(loc="upper left", bbox_to_anchor=(0.08, 1.0))
head(ax, "(a) Subcell J-V characteristics")
ax = axs[1]
# below the last computed point the top cell is reverse-biased: the series current stays at its J_sc
# (R_sh = 1e5 Ohm cm2 adds only (V0 - V)/R_sh), so the 2T curve is continued to 0 V
t2v = [r["V"] for r in t2]; t2j = [r["J_mAcm2"] for r in t2]
k0 = int(np.argmin(t2v)); v0, j0 = t2v[k0], t2j[k0]
ext_v = list(np.linspace(0, v0, 50)); ext_j = [j0 + 1000 * (v0 - v) / 1e5 for v in ext_v]
order = np.argsort(t2v)
ax.plot(ext_v[:-1] + [t2v[k] for k in order], ext_j[:-1] + [t2j[k] for k in order], color=C["2T"], lw=1.6,
        label="2T output")
vv = [r["V"] for r in vm]
ax.plot(vv, [r["J_out_mAcm2"] for r in vm], color=C["VMM"], lw=1.6, label=fr"VMM output ($n$ = {NOPT:.3f})")
ax.plot(vv, [r["J_top_mAcm2"] for r in vm], color=C["VMM"], ls=(0, (5, 2)), lw=1.1, label="VMM perovskite branch")
ax.plot(vv, [r["J_str_mAcm2"] for r in vm], color=C["VMM"], ls=(0, (1.5, 1.5)), lw=1.1, label="VMM Si string")
ax.plot([I["2T"]["Vout"], I["VMM"]["Vout"]], [0.1 * I["2T"]["J"], 0.1 * I["VMM"]["Jout"]], "o", color="black",
        mfc="white", ms=5.5, mew=1.1)
ax.set(xlim=(0, 2.1), ylim=(0, 49), xlabel="Output voltage (V)", ylabel=r"Current density (mA cm$^{-2}$)")
ax.legend(loc="upper right")
head(ax, "(b) 2T and VMM output characteristics")
save(fig, "Fig3_reference_JV")

# ---------------------------------------------------------------- Fig: loss decomposition
cfgs = ["2T", "VMM", "VMM-n2", "3T-L1", "3T-L4", "4T"]
labels = ["2T series", r"VMM, $n^{*}$", r"VMM, $n$ = 2", "3T-L1", "3T-L4", "4T"]
cons = [P0 - 100 * I[c]["P"] / PIN for c in cfgs]
resl = [100 * (I[c]["P"] - ref["Nominal"][c]["P"]) / PIN for c in cfgs]
fig, ax = plt.subplots(figsize=(W2 * 0.78, 6.0 * CM), constrained_layout=True)
y = np.arange(len(cfgs))[::-1]
ax.barh(y, cons, height=0.6, color="#9ec5f4", edgecolor="black", linewidth=0.7, label="Constraint loss ($R$ = 0)")
ax.barh(y, resl, left=cons, height=0.6, color="#1f3fd8", edgecolor="black", linewidth=0.7,
        label=r"Resistive loss ($R_x$ = $R_\ell$ = 0.5 $\Omega$ cm$^2$)")
for yi, a, b in zip(y, cons, resl):
    ax.text(a + b + 0.03, yi, f"{a + b:.2f}", va="center", fontsize=8)
ax.set_yticks(y, labels)
ax.yaxis.set_minor_locator(NullLocator())
ax.tick_params(axis="y", which="both", right=False, left=False)
ax.set_xlabel(r"Connection loss $\Delta\eta_{\rm conn}$ (pp)")
ax.set_xlim(0, 2.2); ax.set_ylim(-0.6, len(cfgs) - 0.4)
ax.legend(loc="lower right")
save(fig, "Fig4_loss_decomposition")

# ---------------------------------------------------------------- Fig: resistance sweeps
rs = read("resistance_sweeps_reference")
fig, axs = plt.subplots(1, 2, figsize=(W2, 6.2 * CM), sharey=True, constrained_layout=True)
for ax, case, txt, pl in zip(axs, ["A", "B"], [r"Interconnection element $R_x$", r"Every terminal lead $R_\ell$"], "ab"):
    sub = [r for r in rs if r["Case"] == case]
    R = [r["R_ohm_cm2"] for r in sub]
    for c in ["2T", "VMM", "3T-L1", "3T-L4", "4T"]:
        line(ax, R, [r["PCE_" + c] for r in sub], c, ms=3.6)
    ax.set_xlabel(r"Area-specific resistance ($\Omega$ cm$^2$)")
    ax.set_xlim(-0.1, 5.1); ax.set_ylim(28, 38.6)
    if pl == "a":
        head(ax, r"(a) Interconnection element $\mathbf{R_x}$ only")
    else:
        head(ax, r"(b) Every terminal lead $\mathbf{R_\ell}$")
axs[0].set_ylabel("Tandem PCE (%)")
axs[1].legend(loc="lower left", ncol=1)
save(fig, "Fig5_resistance_sweeps")

# ---------------------------------------------------------------- Fig: 3T shared contact
mp = read("3T_Rz_vs_mismatch_map")
fig, axs = plt.subplots(1, 2, figsize=(W2, 6.2 * CM), constrained_layout=True)
ax = axs[0]
sel = [0.85, 0.90, 0.95, 1.00, 1.05, 1.10]
cols = ["#0b1f3d", "#163f7a", "#256abf", "#3f86e0", "#6da7ec", "#9ec5f4"]
mks = ["o", "s", "^", "v", "D", "<"]
for f, col, mk in zip(sel, cols, mks):
    sub = [r for r in mp if abs(r["ftop"] - f) < 1e-6]
    ax.plot([r["Rz_ohm_cm2"] for r in sub], [r["Loss_pp"] for r in sub], color=col, marker=mk, mfc=col, mec=col,
            ms=3.6, label=fr"$f_{{\rm QIBC}}$ = {f:.2f}")
ax.set(xlabel=r"Shared-contact resistance $R_Z$ ($\Omega$ cm$^2$)", ylabel="Loss relative to lossless 4T (pp)",
       xlim=(-0.5, 20.5), ylim=(-0.008, 0.23))
ax.legend(loc="upper left", bbox_to_anchor=(0.06, 0.99))
head(ax, "(a) Loss versus shared-contact resistance")
ax = axs[1]
x = np.array([r["Rz_DeltaJ2_pp"] for r in mp]); yv = np.array([r["Loss_pp"] for r in mp])
lim = 0.24
ax.plot([0, lim], [0, lim], color="black", lw=1.0, ls=(0, (5, 2)), label=r"$R_Z\,\Delta J_{\rm mp}^2/P_{\rm in}$")
ax.plot(x, yv, "o", color=C["3T"], mfc="white", mec=C["3T"], ms=4.2, mew=0.9, label="Full re-optimization")
ax.set(xlabel=r"Estimate $R_Z\,\Delta J_{\rm mp}^2/P_{\rm in}$ (pp)", ylabel="Computed loss (pp)", xlim=(0, lim), ylim=(0, lim))
ax.legend(loc="upper left", bbox_to_anchor=(0.06, 0.99))
head(ax, "(b) Computed versus estimated loss")
save(fig, "Fig6_3T_shared_contact")

# ---------------------------------------------------------------- Fig: ratio sweeps
vr = read("vmm_ratio_sweep_reference"); tr = read("3T_L4_ratio_sweep_reference")
fig, ax = plt.subplots(figsize=(W2 * 0.72, 6.4 * CM), constrained_layout=True)
nn = [r["n"] for r in vr]
ax.axhline(P0, color="black", lw=0.9, ls=(0, (5, 2)))
ax.text(1.415, P0 + 0.08, "Lossless 4T", fontsize=8)
ax.axvline(NOPT, color="0.35", lw=0.8, ls=(0, (1.5, 1.5)))
ax.text(NOPT + 0.012, 33.25, fr"$n^{{*}}$ = {NOPT:.3f}", fontsize=8)
ax.plot(nn, [r["PCE_R0"] for r in vr], color=C["VMM"], lw=1.5, label=r"VMM, $R$ = 0")
ax.plot(nn, [r["PCE_Rstr_1ohmcm2"] for r in vr], color=C["VMM"], lw=1.3, ls=(0, (5, 2)), label=r"VMM, $R_{\rm str}$ = 1 $\Omega$ cm$^2$")
ax.plot(nn, [r["PCE_Rtop_1ohmcm2"] for r in vr], color=C["VMM"], lw=1.3, ls=(0, (1.5, 1.5)), label=r"VMM, $R_{\rm QIBC}$ = 1 $\Omega$ cm$^2$")
mm = [r["m"] for r in tr][::2]; pm = [r["PCE_R0"] for r in tr][::2]
ax.plot(mm, pm, "^", color=C["3T"], mfc="white", mec=C["3T"], ms=4.5, mew=0.9, label=r"3T-L4 versus $m$, $R$ = 0")
ax.set(xlabel=r"Voltage-matching ratio $n$ (VMM) or $m$ (3T-L4)", ylabel="Tandem PCE (%)", xlim=(1.4, 2.4), ylim=(33, 38.4))
ax.legend(loc="lower right")
save(fig, "Fig7_ratio_sweeps")

# ---------------------------------------------------------------- Fig: sensitivity (2 x 2)
det = read("sensitivity_detail")
apes = {r["alpha"]: r["APE_350_1050nm_eV"] for r in read("spectral_tilt_APE")}
sweeps = [("spectrum", "Average photon energy (eV)", lambda v: apes[v]),
          ("temperature", "Cell temperature (°C)", lambda v: v),
          ("mismatch", r"Top-cell generation factor $f_{\rm QIBC}$", lambda v: v),
          ("voltage", r"Top-cell $V_{\rm oc}$ (V)", None)]


def sens_fig(rstate, metric, cfgset, ylabel, name, ycap=None):
    fig, axs = plt.subplots(2, 2, figsize=(W2, 11.8 * CM), constrained_layout=True)
    for ax, (sw, xl, xf), pl in zip(axs.flat, sweeps, "abcd"):
        offs, top = [], 0
        for c in cfgset:
            sub = sorted([r for r in det if r["Sweep"] == sw and r["Rstate"] == rstate and r["Config"] == c],
                         key=lambda r: r["Value"])
            xs = [r["Voc_top"] if xf is None else xf(r["Value"]) for r in sub]
            ys = [r[metric] for r in sub]
            line(ax, xs, ys, c, ms=3.3)
            if ycap is not None and max(ys) > ycap:
                offs.append(STYLE[c][0]); top = max(top, max(ys))
        ax.set_xlabel(xl); ax.set_ylabel(ylabel)
        if ycap is not None:
            ax.set_ylim(-0.15, ycap)
            if offs:
                nx, ny, ha = {"voltage": (0.03, 0.42, "left"), "temperature": (0.03, 0.95, "left")}.get(sw, (0.97, 0.62, "right"))
                note(ax, f"Off scale: {', '.join(offs)}\n(maximum {top:.1f} pp)", nx, ny, ha=ha, size=7.5)
        head(ax, f"({pl}) " + {"spectrum": "Spectral variation", "temperature": "Cell temperature",
                               "mismatch": "Subcell current mismatch", "voltage": "Top-cell voltage loss"}[sw])
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=min(6, len(l)), fontsize=7.5)
    save(fig, name)


sens_fig("R0", "Loss_vs_4T0_pp", ["2T", "VMM", "VMM-n2", "3T-L2", "3T-L2b", "3T-L3"],
         "Loss relative to lossless 4T (pp)", "Fig8_sensitivity_constraint_loss", ycap=4.0)
sens_fig("Rnom", "PCE", ["2T", "VMM", "3T-L1", "3T-L2", "3T-L4", "4T"], "Tandem PCE (%)", "Fig9_sensitivity_nominal_R")

# ---------------------------------------------------------------- Fig: design windows (2 x 2)
def window(ax, rows, xk, yk, zk, xl, yl, vmax=2.0):
    xs = sorted({r[xk] for r in rows}); ys = sorted({r[yk] for r in rows})
    Z = np.full((len(ys), len(xs)), np.nan)
    for r in rows:
        Z[ys.index(r[yk]), xs.index(r[xk])] = r[zk]
    levels = np.linspace(0, vmax, 21)
    im = ax.contourf(xs, ys, np.clip(Z, 0, vmax), levels=levels, cmap=SEQ_CMAP, extend="max")
    cs = ax.contour(xs, ys, Z, levels=[0.1, 0.25, 0.5, 1.0], colors="black", linewidths=0.7)
    ax.clabel(cs, fmt=lambda v: f"{v:g}", fontsize=7, inline_spacing=2)
    ax.set(xlabel=xl, ylabel=yl)
    return im


fig, axs = plt.subplots(2, 2, figsize=(W2, 11.6 * CM), constrained_layout=True)
im = window(axs[0, 0], read("VMM_n_vs_mismatch_map"), "n", "ftop", "Loss_pp", r"VMM ratio $n$", r"Top-cell generation factor $f_{\rm QIBC}$")
window(axs[0, 1], read("VMM_n_vs_temperature_map"), "n", "T_C", "Loss_pp", r"VMM ratio $n$", "Cell temperature (°C)")
window(axs[1, 0], [dict(r, ape=r["APE_eV"]) for r in read("2T_alpha_vs_mismatch_map")], "ape", "ftop", "Loss_pp",
       "Average photon energy (eV)", r"Top-cell generation factor $f_{\rm QIBC}$")
window(axs[1, 1], read("3T_Rz_vs_mismatch_map"), "Rz_ohm_cm2", "ftop", "Loss_pp",
       r"Shared-contact resistance $R_Z$ ($\Omega$ cm$^2$)", r"Top-cell generation factor $f_{\rm QIBC}$")
for ax, pl, lab in zip(axs.flat, "abcd", ["VMM: ratio vs. mismatch", "VMM: ratio vs. temperature",
                                         "2T: spectrum vs. mismatch", r"3T-L1: $\mathbf{R_Z}$ vs. mismatch"]):
    head(ax, f"({pl}) {lab}")
cb = fig.colorbar(im, ax=axs, shrink=0.75, pad=0.02, ticks=[0, 0.5, 1.0, 1.5, 2.0])
cb.set_label("Loss relative to lossless 4T (pp)")
cb.outline.set_linewidth(0.9)
cb.ax.tick_params(direction="in", which="both")
save(fig, "Fig11_design_windows")

# ---------------------------------------------------------------- Fig: factorial robustness (shortfall)
fac = read("factorial_robustness")
show = ["2T", "VMM", "VMM-n2", "3T-L1", "3T-L2", "3T-L2b", "3T-L3", "3T-L4", "4T"]
fig, axs = plt.subplots(1, 2, figsize=(W2, 7.0 * CM), sharey=True, constrained_layout=True)
rng = np.random.default_rng(3)
for ax, st, txt, pl in zip(axs, ["R0", "Rnom"], [r"$R$ = 0", r"$R_x$ = $R_\ell$ = 0.5 $\Omega$ cm$^2$"], "ab"):
    sub = [r for r in fac if r["Rstate"] == st]
    best = np.array([max(r["PCE_" + c] for c in show) for r in sub])
    for i, c in enumerate(show):
        v = best - np.array([r["PCE_" + c] for r in sub])
        lab, col, ls, mk = STYLE[c]
        ax.plot(rng.uniform(-0.2, 0.2, len(v)) + i, v, mk, color=col, mfc="none", mec=col, ms=3.6, mew=0.8)
        ax.plot([i - 0.32, i + 0.32], [v.mean()] * 2, color="black", lw=1.5)
    ax.set_yscale("symlog", linthresh=0.1, linscale=0.6)
    ax.set_xticks(range(len(show)), [STYLE[c][0] for c in show], rotation=90)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(-0.6, len(show) - 0.4); ax.set_ylim(-0.02, 15)
    ax.yaxis.set_major_locator(FixedLocator([0, 0.1, 0.3, 1, 3, 10]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels(["0", "0.1", "0.3", "1", "3", "10"])
    head(ax, f"({pl}) " + {"a": r"Ideal state, $\mathbf{R}$ = 0",
                           "b": r"Nominal state, $\mathbf{R_x = R_\ell}$ = 0.5 $\mathbf{\Omega}$ cm$^\mathbf{2}$"}[pl])
axs[0].set_ylabel("Shortfall from best configuration (pp)")
save(fig, "Fig10_factorial_robustness")
print("figures written to", FIG)
