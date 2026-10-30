"""Robustness metrics and design-window limits from the CSV outputs of wiring_analysis.wls."""
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parent / "Results"


def read(name):
    with open(RES / f"{name}.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    out = []
    for r in rows:
        d = {}
        for k, v in r.items():
            try:
                d[k] = float(v)
            except (TypeError, ValueError):
                d[k] = v
        out.append(d)
    return out


CFG = ["2T", "VMM", "VMM-n2", "3T-L1", "3T-L2", "3T-L2b", "3T-L3", "3T-L4", "4T"]
summary = {}

# ---- factorial robustness metrics
fac = read("factorial_robustness")
rob = {}
for st in ["R0", "Rnom"]:
    sub = [r for r in fac if r["Rstate"] == st]
    best = np.array([max(r["PCE_" + c] for c in CFG) for r in sub])
    ref4 = np.array([r["PCE_4T_lossless_ref"] for r in sub])
    rob[st] = {}
    for c in CFG:
        v = np.array([r["PCE_" + c] for r in sub])
        regret = best - v
        rob[st][c] = {
            "mean": v.mean(), "min": v.min(), "max": v.max(), "std": v.std(),
            "mean_loss_vs_4T0": (ref4 - v).mean(), "max_loss_vs_4T0": (ref4 - v).max(),
            "mean_regret": regret.mean(), "max_regret": regret.max(),
            "share_within_0p1pp_of_best": float((regret <= 0.1).mean()),
            "share_best": float((regret <= 1e-6).mean()),
        }
summary["factorial"] = rob
summary["factorial_n_conditions"] = len([r for r in fac if r["Rstate"] == "R0"])

# ---- one-at-a-time ranges: max loss vs lossless 4T in each sweep
det = read("sensitivity_detail")
oat = defaultdict(dict)
for sw in ["spectrum", "temperature", "mismatch", "voltage"]:
    for st in ["R0", "Rnom"]:
        for c in CFG:
            sub = [r for r in det if r["Sweep"] == sw and r["Rstate"] == st and r["Config"] == c]
            if sub:
                oat[f"{sw}|{st}"][c] = {"max_loss": max(r["Loss_vs_4T0_pp"] for r in sub),
                                         "min_pce": min(r["PCE"] for r in sub),
                                         "max_pce": max(r["PCE"] for r in sub)}
summary["one_at_a_time"] = oat
summary["voltage_sweep_voc_top"] = sorted({(r["Value"], r["Voc_top"]) for r in det if r["Sweep"] == "voltage"})
summary["temperature_voc"] = sorted({(r["Value"], r["Voc_top"], r["Voc_bot"]) for r in det if r["Sweep"] == "temperature"})
summary["spectrum_currents"] = sorted({(r["alpha"], r["APE_eV"], r["Jg_top_mAcm2"], r["Jg_bot_mAcm2"])
                                       for r in read("sensitivity_one_at_a_time") if r["Sweep"] == "spectrum"})


def window_1d(xs, losses, thr):
    ok = [x for x, l in zip(xs, losses) if l <= thr]
    return (min(ok), max(ok)) if ok else None


# ---- VMM n windows
vm = read("VMM_n_vs_mismatch_map")
vmT = read("VMM_n_vs_temperature_map")
win = {}
for thr in (0.1, 0.25, 0.5):
    # n window satisfied at every ftop within 0.95-1.05 and every T 15-75 C
    ns = sorted({r["n"] for r in vm})
    okf = [n for n in ns if all(r["Loss_pp"] <= thr for r in vm if r["n"] == n and 0.949 <= r["ftop"] <= 1.051)]
    okT = [n for n in ns if all(r["Loss_pp"] <= thr for r in vmT if r["n"] == n)]
    okT45 = [n for n in ns if all(r["Loss_pp"] <= thr for r in vmT if r["n"] == n and r["T_C"] <= 45.1)]
    win[f"VMM_n_window_ftop0.95-1.05_{thr}pp"] = (min(okf), max(okf)) if okf else None
    win[f"VMM_n_window_T15-75_{thr}pp"] = (min(okT), max(okT)) if okT else None
    win[f"VMM_n_window_T15-45_{thr}pp"] = (min(okT45), max(okT45)) if okT45 else None
    ref_rows = sorted([r for r in vm if abs(r["ftop"] - 1) < 1e-9], key=lambda r: r["n"])
    win[f"VMM_n_window_reference_{thr}pp"] = window_1d([r["n"] for r in ref_rows], [r["Loss_pp"] for r in ref_rows], thr)
# best fixed n over the temperature range (min of max loss)
ns = sorted({r["n"] for r in vmT})
worst = {n: max(r["Loss_pp"] for r in vmT if r["n"] == n) for n in ns}
nbest = min(worst, key=worst.get)
win["VMM_minimax_n_T15-75"] = (nbest, worst[nbest])
win["VMM_n_opt_by_T"] = [(t, min((r for r in vmT if r["T_C"] == t), key=lambda r: r["Loss_pp"])["n"])
                         for t in sorted({r["T_C"] for r in vmT})]

# ---- 2T windows
t2 = read("2T_alpha_vs_mismatch_map")
ref2 = sorted([r for r in t2 if abs(r["alpha"]) < 1e-9], key=lambda r: r["ftop"])
for thr in (0.1, 0.25, 0.5, 1.0):
    win[f"2T_ftop_window_AM15G_{thr}pp"] = window_1d([r["ftop"] for r in ref2], [r["Loss_pp"] for r in ref2], thr)
    f1 = sorted([r for r in t2 if abs(r["ftop"] - 1) < 1e-9], key=lambda r: r["APE_eV"])
    win[f"2T_APE_window_ftop1_{thr}pp"] = window_1d([r["APE_eV"] for r in f1], [r["Loss_pp"] for r in f1], thr)
win["2T_best_ftop_by_APE"] = [(a, min((r for r in t2 if r["APE_eV"] == a), key=lambda r: r["Loss_pp"])["ftop"])
                              for a in sorted({r["APE_eV"] for r in t2})]

# ---- 3T shared contact
mp = read("3T_Rz_vs_mismatch_map")
win["3T_max_loss_Rz20_any_ftop"] = max(r["Loss_pp"] for r in mp if r["Rz_ohm_cm2"] == 20)
win["3T_max_loss_Rz5_any_ftop"] = max(r["Loss_pp"] for r in mp if r["Rz_ohm_cm2"] == 5)
win["3T_DeltaJ_range"] = (min(r["DeltaJ_mAcm2"] for r in mp), max(r["DeltaJ_mAcm2"] for r in mp))
summary["windows"] = win

# ---- resistance slopes (pp per ohm cm^2, from 0 -> 1 ohm cm^2)
rs = read("resistance_sweeps_reference")
slopes = {}
for case in "AB":
    r0 = next(r for r in rs if r["Case"] == case and r["R_ohm_cm2"] == 0)
    r1 = next(r for r in rs if r["Case"] == case and r["R_ohm_cm2"] == 1)
    slopes[case] = {c: r0["PCE_" + c] - r1["PCE_" + c] for c in ["2T", "VMM", "3T-L1", "3T-L4", "4T"]}
summary["resistance_loss_at_1ohmcm2"] = slopes

summary["reference"] = json.load(open(RES / "reference_operating_points.json"))
summary["vmm_branch"] = read("vmm_branch_resistance")
summary["verification"] = read("verification")
json.dump(summary, open(RES / "chapter_summary.json", "w"), indent=1, default=float)


def show(d, ind=0):
    for k, v in d.items():
        if isinstance(v, dict):
            print(" " * ind + str(k)); show(v, ind + 2)
        else:
            print(" " * ind + f"{k}: {v}")


show({"factorial": rob, "windows": win, "slopes": slopes})
for k, v in oat.items():
    print(k, {c: round(x["max_loss"], 3) for c, x in v.items()})
print("voc_top vs kNR", summary["voltage_sweep_voc_top"])
print("T voc", summary["temperature_voc"])
print("spectrum", summary["spectrum_currents"])
