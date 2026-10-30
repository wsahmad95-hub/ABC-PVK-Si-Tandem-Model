(* ::Package:: *)
(* Wiring analysis: interconnection losses and robustness of 2T, VMM, 3T and 4T
   all-back-contact perovskite/Si tandems.
   Uses the subcell equations, parameters and optics of ABC-PVK-Si_Tandem_Model.nb
   (state saved by run_model.wls in model_state.mx). Nothing in the subcell
   physics is changed; only the circuit layer and the operating conditions are varied.

   Units: internal current density A/m^2, voltage V, resistance input in ohm cm^2
   (converted with 1e-4 to ohm m^2), power density W/m^2.                          *)

work = DirectoryName[$InputFileName] <> "/";
Get[work <> "model_state.mx"];
out = work <> "Results/";
If[! DirectoryQ[out], CreateDirectory[out]];
$HistoryLength = 0;
tic = AbsoluteTime[];
log[x__] := (Print[Round[AbsoluteTime[] - tic, 0.1], " s | ", x]; );

Tref = 298.15; Pin = N[PowerofSun];
c2 = 10.^-4;  (* ohm cm^2 -> ohm m^2 *)

(* ---------------------------------------------------------------- optics *)
safe[f_InterpolatingFunction] := Module[{d = f["Domain"][[1]]},
   Function[x, If[d[[1]] <= x <= d[[2]], Max[0., f[x]], 0.]]];
eqePe = safe[EQEPerovskite]; eqeSi = safe[EQESi]; absTop = safe[TransmittedtoSi];
gamma0 = Interpolation[makeFlux[AM1p5global]];
gdom = gamma0["Domain"][[1]];
dE = 0.0005;
egrid = Range[Ceiling[gdom[[1]], dE], 4.4, dE];  (* full spectrum so the tilt keeps total irradiance *)
gam0vals = If[gdom[[1]] <= # <= gdom[[2]], gamma0[#], 0.] & /@ egrid;
P0E = Total[gam0vals egrid] dE;           (* irradiance from the photon-flux table *)
spectrum[alpha_] := Module[{w = (egrid/2.)^alpha, s},
   s = gam0vals w; s P0E/(Total[s egrid] dE)];  (* tilt, renormalised to the same irradiance *)
apeBand = Boole[1.1808 <= # <= 3.5424] & /@ egrid;  (* 350-1050 nm *)
ape[alpha_] := With[{s = spectrum[alpha]}, Total[s egrid apeBand]/Total[s apeBand]];

gens[alpha_, Tk_, ftop_ : 1.] := Module[{s = spectrum[alpha], dEp, dEs, top, bot, e1, e2},
   dEp = PerovskiteBandgap[Tk] - PerovskiteBandgap[Tref];
   dEs = SiBandgap[Tk] - SiBandgap[Tref];
   e1 = Boole[1.2 <= # <= 4.4] & /@ egrid; e2 = Boole[0.8 <= # <= 3.5] & /@ egrid;
   top = Total[(eqePe[# - dEp] & /@ egrid) s e1] dE;
   bot = Total[(Max[0., eqeSi[# - dEs] - absTop[# - dEp]] & /@ egrid) s e2] dE;
   {ftop top, bot}];

(* ------------------------------------------------------- temperature state *)
setT[Tk_] := (T = N[Tk];
   JRad0Si = q 2 Pi/(c^2 h2^3) NIntegrate[En^2/(Exp[En/(k T)] - 1), {En, SiBandgap[T], 4.4}];
   JRad0Perovskite = q 2 Pi/(c^2 h2^3) NIntegrate[En^2/(Exp[En/(k T)] - 1), {En, PerovskiteBandgap[T], 4.4}];);

(* ------------------------------------------------------------ subcell J-V *)
buildCurve[type_String, g_?NumericQ, kNR_ : 1., dv_ : 0.001] := Module[
   {res, voc, vmax, grid, jat, data, fwd, z, y, egap},
   res = If[type === "Pe",
     Function[{x, j}, j - g + JRPerovskite[x, j, 2, RsPerovskite] +
        kNR JRnonPerovskite[x, j, RsPerovskite] + Resistance[x, j, RsPerovskite, RshPerovskite]],
     Function[{x, j}, j - g + JRSi[x, j, 1, RsSi] + JRAuger[x, j, RsSi] +
        JRnonSi[x, j, RsSi] + Resistance[x, j, RsSi, RshSi]]];
   egap = If[type === "Pe", PerovskiteBandgap[T], SiBandgap[T]];
   voc = z /. FindRoot[res[z, 0.] == 0, {z, 0.2, egap}, Method -> "Brent", AccuracyGoal -> 12];
   vmax = voc + If[type === "Pe", 0.2, 0.3];
   jat[x_] := Module[{lo = -N[g], hi = 1.01 N[g] + 1., cnt = 0},
     While[res[x, lo] > 0 && cnt < 80, lo = 2 lo; cnt++];
     y /. FindRoot[res[x, y] == 0, {y, lo, hi}, Method -> "Brent", AccuracyGoal -> 10, PrecisionGoal -> 10]];
   grid = Union[Range[0., vmax, dv], {voc}];
   data = Quiet[{#, If[# == voc, 0., jat[#]]} & /@ grid];
   fwd = SortBy[Reverse /@ Select[data, #[[2]] >= 0 &], First];
   <|"Type" -> type, "g" -> g, "Voc" -> voc, "Jsc" -> data[[1, 2]], "Vmax" -> Last[grid], "Data" -> data,
     "J" -> Interpolation[data, InterpolationOrder -> 3],
     "V" -> Interpolation[fwd, InterpolationOrder -> 3]|>];

curvesAt[alpha_, Tk_, ftop_ : 1., kNR_ : 1., dv_ : 0.001] := Module[{gp, gs},
   setT[Tk]; {gp, gs} = gens[alpha, Tk, ftop];
   {buildCurve["Pe", gp, kNR, dv], buildCurve["Si", gs, 1., dv]}];

jf[c_, v_?NumericQ] := c["J"][Clip[v, {0., c["Vmax"]}]];

(* ------------------------------------------------------------ maximizers *)
(* grid search followed by golden-section refinement in the bracketing cell *)
numOr[v_] := If[NumericQ[v], N[v], -10.^9];
max1D[f_, a_, b_, n_ : 300] := Module[{xs, vals, i, lo, hi, gr = (Sqrt[5.] - 1)/2, x1, x2, f1, f2, ff, best},
   ff = Function[arg, numOr[Quiet[f[arg]]]];
   xs = Subdivide[N[a], N[b], n]; vals = ff /@ xs;
   i = First@Ordering[vals, -1];
   lo = xs[[Max[1, i - 1]]]; hi = xs[[Min[n + 1, i + 1]]];
   x1 = hi - gr (hi - lo); x2 = lo + gr (hi - lo); f1 = ff[x1]; f2 = ff[x2];
   Do[If[f1 >= f2, hi = x2; x2 = x1; f2 = f1; x1 = hi - gr (hi - lo); f1 = ff[x1],
      lo = x1; x1 = x2; f1 = f2; x2 = lo + gr (hi - lo); f2 = ff[x2]];
    If[hi - lo < 10^-11, Break[]], {60}];
   best = If[f1 >= f2, {x1, f1}, {x2, f2}];
   If[best[[2]] >= vals[[i]], best, {xs[[i]], vals[[i]]}]];

(* ------------------------------------------------------- configurations
   All resistances in ohm cm^2 of the common aperture.                     *)

(* 4T: each port has its own series path resistance Rport *)
port1[cur_, R_] := Module[{r = c2 R, s},
   s = max1D[Function[v, With[{j = jf[cur, v]}, v j - r j^2]], 0., cur["Voc"]];
   <|"V" -> s[[1]] - r jf[cur, s[[1]]], "J" -> jf[cur, s[[1]]], "P" -> s[[2]]|>];
config4T[pe_, si_, Rp_, Rs_] := Module[{a = port1[pe, Rp], b = port1[si, Rs]},
   <|"P" -> a["P"] + b["P"], "VTop" -> a["V"], "VBot" -> b["V"], "JTop" -> a["J"], "JBot" -> b["J"],
     "Loss_R" -> c2 (Rp a["J"]^2 + Rs b["J"]^2)|>];

(* 2T: one common current; Rser lumps interconnect and both leads *)
config2T[pe_, si_, Rser_] := Module[{r = c2 Rser, jmax, s, j},
   jmax = 0.999999 Min[pe["Jsc"], si["Jsc"]];
   s = max1D[Function[x, x (pe["V"][x] + si["V"][x] - r x)], 0., jmax, 400];
   j = s[[1]];
   <|"P" -> s[[2]], "J" -> j, "VPe" -> pe["V"][j], "VSi" -> si["V"][j], "Vout" -> pe["V"][j] + si["V"][j] - r j,
     "Loss_R" -> r j^2, "LimitingCell" -> If[pe["Jsc"] < si["Jsc"], "Pe", "Si"]|>];

(* VMM: one QIBC branch in parallel with a string of n Si segments sharing the aperture.
   Segment current = J_Si/n, string voltage = n V_Si. Rtop, Rstr: branch resistances,
   Rbus: common bus/tie carrying the total current. *)
configVMM[pe_, si_, n_, Rtop_, Rstr_, Rbus_] := Module[{top, str, vlo, vhi, jt, s, vb, rb = c2 Rbus},
   top = Interpolation[{#[[1]] - c2 Rtop #[[2]], #[[2]]} & /@ pe["Data"], InterpolationOrder -> 3];
   str = Interpolation[{n #[[1]] - c2 Rstr #[[2]]/n, #[[2]]/n} & /@ si["Data"], InterpolationOrder -> 3];
   vlo = Max[top["Domain"][[1, 1]], str["Domain"][[1, 1]], 0.];
   vhi = Min[top["Domain"][[1, 2]], str["Domain"][[1, 2]]];
   jt[v_?NumericQ] := top[v] + str[v];
   s = max1D[Function[v, (v - rb jt[v]) jt[v]], vlo, vhi, 400];
   vb = s[[1]];
   <|"P" -> s[[2]], "Vbus" -> vb, "Vout" -> vb - rb jt[vb], "JTop" -> top[vb], "JStr" -> str[vb],
     "Jout" -> jt[vb], "VSiCell" -> (vb + c2 Rstr str[vb])/n,
     "Loss_R" -> c2 (Rtop top[vb]^2 + Rstr str[vb]^2) + rb jt[vb]^2,
     "AtDomainEdge" -> Boole[vb > vhi - 0.002 || vb < vlo + 0.002]|>];

(* 3T common-Z s-type (perovskite negative joined to Si positive at node N, Z lead from N).
   Port equations:  V_TZ = vp - Rt Jp - Rz (Jp - Js),   V_ZR = vs - Rr Js + Rz (Jp - Js)
   P = V_TZ Jp + V_ZR Js = vp Jp + vs Js - Rt Jp^2 - Rr Js^2 - Rz (Jp - Js)^2              *)
ports3[pe_, si_, Rt_, Rr_, Rz_, vp_, vs_] := Module[{jp = jf[pe, vp], js = jf[si, vs]},
   {vp - c2 Rt jp - c2 Rz (jp - js), vs - c2 Rr js + c2 Rz (jp - js), jp, js}];
(* each output port feeds a converter input that cannot return power (blocking diode):
   delivered port power is clipped at zero. Only matters for the fixed-voltage loadings L2/L3. *)
pack3[pe_, si_, Rt_, Rr_, Rz_, vp_, vs_] := Module[{p = ports3[pe, si, Rt, Rr, Rz, vp, vs]},
   <|"P" -> Max[0., p[[1]] p[[3]]] + Max[0., p[[2]] p[[4]]], "PunClipped" -> p[[1]] p[[3]] + p[[2]] p[[4]],
     "VTZ" -> p[[1]], "VZR" -> p[[2]], "JTop" -> p[[3]], "JBot" -> p[[4]],
     "JZ" -> p[[4]] - p[[3]], "vp" -> vp, "vs" -> vs,
     "Loss_R" -> c2 (Rt p[[3]]^2 + Rr p[[4]]^2 + Rz (p[[3]] - p[[4]])^2)|>];

(* L1: two independent MPP trackers (joint maximum, coordinate ascent) *)
config3TL1[pe_, si_, Rt_, Rr_, Rz_] := Module[{vp, vs, pw, it, old = -1., new},
   pw[x_, y_] := Module[{p = ports3[pe, si, Rt, Rr, Rz, x, y]}, p[[1]] p[[3]] + p[[2]] p[[4]]];
   vp = port1[pe, 0.]["V"]; vs = port1[si, 0.]["V"];
   Do[
    vp = max1D[Function[x, pw[x, vs]], 0., pe["Voc"], 200][[1]];
    {vs, new} = max1D[Function[y, pw[vp, y]], 0., si["Voc"], 200];
    If[Abs[new - old] < 10^-9, Break[]]; old = new, {it, 12}];
   pack3[pe, si, Rt, Rr, Rz, vp, vs]];

(* inner solve: given vp and a target relation, find vs (monotone in vs) *)
(* bisection: fun is monotone increasing in vs on [0, Vmax] *)
solveVs[si_, fun_] := Module[{lo = 0., hi = N[si["Vmax"]], flo, fhi, mid, fm},
   flo = fun[lo]; fhi = fun[hi];
   If[! (NumericQ[flo] && NumericQ[fhi] && flo <= 0 <= fhi), Return[$Failed]];
   Do[mid = (lo + hi)/2; fm = fun[mid];
    If[fm > 0, hi = mid, lo = mid]; If[hi - lo < 10^-12, Break[]], {60}];
   (lo + hi)/2];

(* L2: TZ port tracked, ZR port held at a fixed external voltage VZRfix *)
config3TL2[pe_, si_, Rt_, Rr_, Rz_, VZRfix_] := Module[{vsOf, pw, s},
   vsOf[vp_] := Module[{jp = jf[pe, vp]},
     solveVs[si, Function[y, y - c2 Rr jf[si, y] + c2 Rz (jp - jf[si, y]) - VZRfix]]];
   pw[vp_] := Module[{vs = vsOf[vp]}, If[vs === $Failed, -10.^9, pack3[pe, si, Rt, Rr, Rz, vp, vs]["P"]]];
   s = max1D[pw, 0., pe["Voc"], 200];
   pack3[pe, si, Rt, Rr, Rz, s[[1]], vsOf[s[[1]]]]];

(* L2b: ZR port tracked, TZ (perovskite) port held at a fixed external voltage VTZfix *)
config3TL2b[pe_, si_, Rt_, Rr_, Rz_, VTZfix_] := Module[{vpOf, pw, s},
   vpOf[vs_] := Module[{js = jf[si, vs]},
     solveVs[pe, Function[x, x - c2 Rt jf[pe, x] - c2 Rz (jf[pe, x] - js) - VTZfix]]];
   pw[vs_] := Module[{vp = vpOf[vs]}, If[vp === $Failed, -10.^9, pack3[pe, si, Rt, Rr, Rz, vp, vs]["P"]]];
   s = max1D[pw, 0., si["Voc"], 200];
   pack3[pe, si, Rt, Rr, Rz, vpOf[s[[1]]], s[[1]]]];

(* L3: both port voltages fixed (no tracking) *)
config3TL3[pe_, si_, Rt_, Rr_, Rz_, VTZfix_, VZRfix_] := Module[{sol, x, y, e1, e2},
   e1[a_?NumericQ, b_?NumericQ] := ports3[pe, si, Rt, Rr, Rz, a, b][[1]] - VTZfix;
   e2[a_?NumericQ, b_?NumericQ] := ports3[pe, si, Rt, Rr, Rz, a, b][[2]] - VZRfix;
   sol = Quiet@FindRoot[{e1[x, y] == 0, e2[x, y] == 0}, {{x, VTZfix}, {y, VZRfix}}];
   pack3[pe, si, Rt, Rr, Rz, x /. sol, y /. sol]];

(* L4: single tracker with fixed port-voltage ratio V_TZ = m V_ZR (3T voltage-matched loading) *)
config3TL4[pe_, si_, Rt_, Rr_, Rz_, m_] := Module[{vsOf, pw, s},
   vsOf[vp_] := Module[{jp = jf[pe, vp]},
     solveVs[si, Function[y, Module[{js = jf[si, y]},
        (y - c2 Rr js + c2 Rz (jp - js)) - (vp - c2 Rt jp - c2 Rz (jp - js))/m]]]];
   pw[vp_] := Module[{vs = vsOf[vp]}, If[vs === $Failed, -10.^9, pack3[pe, si, Rt, Rr, Rz, vp, vs]["P"]]];
   s = max1D[pw, 0., pe["Voc"], 200];
   pack3[pe, si, Rt, Rr, Rz, s[[1]], vsOf[s[[1]]]]];

(* ----------------------------------------------- resistance state mapping
   Rx: interconnection element created by the wiring (2T series link P-/S+, VMM parallel tie/bus,
       3T shared Z contact; 4T has none).
   Rl: identical resistance of every external terminal lead (2T: 2 leads, VMM: 2, 3T: 3, 4T: 4).  *)
evalAll[pe_, si_, st_Association] := Module[{rx = st["Rx"], rl = st["Rl"], res = <||>},
   res["4T"] = config4T[pe, si, 2 rl, 2 rl];
   res["2T"] = config2T[pe, si, rx + 2 rl];
   res["VMM"] = configVMM[pe, si, st["n"], 0., 0., rx + 2 rl];
   If[KeyExistsQ[st, "n2"], res["VMM-n2"] = configVMM[pe, si, st["n2"], 0., 0., rx + 2 rl]];
   res["3T-L1"] = config3TL1[pe, si, rl, rl, rx + rl];
   If[KeyExistsQ[st, "VZRfix"],
    res["3T-L2"] = config3TL2[pe, si, rl, rl, rx + rl, st["VZRfix"]];
    res["3T-L2b"] = config3TL2b[pe, si, rl, rl, rx + rl, st["VTZfix"]];
    res["3T-L3"] = config3TL3[pe, si, rl, rl, rx + rl, st["VTZfix"], st["VZRfix"]]];
   res["3T-L4"] = config3TL4[pe, si, rl, rl, rx + rl, st["m"]];
   res["4T0"] = config4T[pe, si, 0., 0.];  (* lossless reference: independent MPPs, no added R *)
   res];

pce[p_] := 100. p/Pin;
cfgOrder = {"2T", "VMM", "VMM-n2", "3T-L1", "3T-L2", "3T-L2b", "3T-L3", "3T-L4", "4T"};

rowOf[tag_, cond_List, gp_, gs_, res_] := Join[{tag}, cond, {0.1 gp, 0.1 gs},
   Table[If[KeyExistsQ[res, c], pce[res[c]["P"]], ""], {c, cfgOrder}], {pce[res["4T0"]["P"]]}];
hdr[condNames_] := Join[{"Sweep"}, condNames, {"Jg_top_mAcm2", "Jg_bot_mAcm2"},
   ("PCE_" <> # & /@ cfgOrder), {"PCE_4T_lossless_ref"}];
csv[name_, rows_] := Export[out <> name <> ".csv", rows, "CSV"];

(* ========================================================= 1. REFERENCE *)
log["Reference curves"];
{peR, siR} = curvesAt[0., Tref];
log["Jg top ", 0.1 peR["g"], "  Jg bot ", 0.1 siR["g"], " (notebook ", 0.1 jgPerovskitetandem, ", ", 0.1 jgSitandem, ")"];

(* optimal VMM ratio and 3T ratio at reference, R = 0 *)
nGrid = Range[1.40, 2.40, 0.01];
vmmN = Table[{nn, pce[configVMM[peR, siR, nn, 0., 0., 0.]["P"]]}, {nn, nGrid}];
nOptCoarse = First[MaximalBy[vmmN, Last]][[1]];
nOpt = max1D[Function[nn, configVMM[peR, siR, nn, 0., 0., 0.]["P"]], nOptCoarse - 0.02, nOptCoarse + 0.02, 20][[1]];
log["VMM optimum n = ", nOpt];
m3 = nOpt;  (* 3T-L4 ratio chosen at the same design point *)

stIdeal = <|"Rx" -> 0., "Rl" -> 0., "n" -> nOpt, "n2" -> 2., "m" -> m3|>;
RNOM = 0.5;
stNom = <|"Rx" -> RNOM, "Rl" -> RNOM, "n" -> nOpt, "n2" -> 2., "m" -> m3|>;

(* fixed-voltage references for L2 / L3 taken from the reference L1 point of each R state *)
Do[With[{s = st},
   ref = config3TL1[peR, siR, s["Rl"], s["Rl"], s["Rx"] + s["Rl"]];
   If[s === stIdeal, stIdeal = Join[stIdeal, <|"VTZfix" -> ref["VTZ"], "VZRfix" -> ref["VZR"]|>],
    stNom = Join[stNom, <|"VTZfix" -> ref["VTZ"], "VZRfix" -> ref["VZR"]|>]]], {st, {stIdeal, stNom}}];
log["Fixed voltages ideal: ", {stIdeal["VTZfix"], stIdeal["VZRfix"]}, " nominal: ", {stNom["VTZfix"], stNom["VZRfix"]}];

refIdeal = evalAll[peR, siR, stIdeal];
refNom = evalAll[peR, siR, stNom];
log["Reference ideal: ", pce[#["P"]] & /@ refIdeal];
log["Reference nominal: ", pce[#["P"]] & /@ refNom];

flat[a_Association] := Normal[a /. x_?NumericQ :> N[x]];
Export[out <> "reference_operating_points.json",
  <|"IncidentPower_W_m2" -> Pin, "T_K" -> Tref, "Jg_top_A_m2" -> peR["g"], "Jg_bot_A_m2" -> siR["g"],
    "Voc_top" -> peR["Voc"], "Voc_bot" -> siR["Voc"], "Jsc_top_A_m2" -> peR["Jsc"], "Jsc_bot_A_m2" -> siR["Jsc"],
    "nOpt" -> nOpt, "m3" -> m3, "RNOM_ohm_cm2" -> RNOM,
    "Ideal" -> (Association[Select[Normal[#], NumericQ[Last[#]] &]] & /@ refIdeal),
    "Nominal" -> (Association[Select[Normal[#], NumericQ[Last[#]] &]] & /@ refNom),
    "FixedV_ideal" -> {stIdeal["VTZfix"], stIdeal["VZRfix"]}, "FixedV_nominal" -> {stNom["VTZfix"], stNom["VZRfix"]}|>, "JSON"];

csv["vmm_ratio_sweep_reference", Prepend[
   Table[{nn, pce[configVMM[peR, siR, nn, 0., 0., 0.]["P"]], pce[configVMM[peR, siR, nn, 0., 0., RNOM + 2 RNOM]["P"]],
     pce[configVMM[peR, siR, nn, 0., 1.0, 0.]["P"]], pce[configVMM[peR, siR, nn, 1.0, 0., 0.]["P"]]}, {nn, nGrid}],
   {"n", "PCE_R0", "PCE_nominal_bus", "PCE_Rstr_1ohmcm2", "PCE_Rtop_1ohmcm2"}]];
csv["3T_L4_ratio_sweep_reference", Prepend[
   Table[{mm, pce[config3TL4[peR, siR, 0., 0., 0., mm]["P"]], pce[config3TL4[peR, siR, RNOM, RNOM, 2 RNOM, mm]["P"]]}, {mm, Range[1.5, 2.3, 0.02]}],
   {"m", "PCE_R0", "PCE_nominal"}]];
log["Ratio sweeps done"];

(* reference J-V curves for figures (common optics) *)
csv["reference_subcell_JV", Prepend[Select[Transpose[{peR["Data"][[All, 1]], 0.1 peR["Data"][[All, 2]]}], #[[2]] >= -5 &], {"V", "J_top_mAcm2"}]];
csv["reference_si_JV", Prepend[Select[Transpose[{siR["Data"][[All, 1]], 0.1 siR["Data"][[All, 2]]}], #[[2]] >= -5 &], {"V", "J_bot_mAcm2"}]];
csv["reference_2T_JV", Prepend[Table[{peR["V"][j] + siR["V"][j], 0.1 j}, {j, Subdivide[0., 0.999999 Min[peR["Jsc"], siR["Jsc"]], 400]}], {"V", "J_mAcm2"}]];
csv["reference_VMM_JV", Prepend[Table[{v, 0.1 jf[peR, v], 0.1 jf[siR, v/nOpt]/nOpt, 0.1 (jf[peR, v] + jf[siR, v/nOpt]/nOpt)}, {v, 0., 1.36, 0.002}],
   {"V", "J_top_mAcm2", "J_str_mAcm2", "J_out_mAcm2"}]];

(* ===================================================== 2. VERIFICATION *)
log["Verification against notebook MakeTandemJVPlots"];
setT[Tref];
verRows = {{"Check", "Notebook_value", "Chapter_code_value", "Difference"}};
nbA = MakeTandemJVPlots["SeriesOptics" -> "SameAsOtherTandems", "n" -> 1.88, "Points" -> 1200];
mine = <|"2T" -> config2T[peR, siR, 0.], "VMM" -> configVMM[peR, siR, 1.88, 0., 0., 0.],
   "3T" -> config3TL1[peR, siR, 0., 0., 0.], "4T" -> config4T[peR, siR, 0., 0.]|>;
Do[AppendTo[verRows, {"R=0 common optics PCE " <> nm, nbA["SummaryData"][[i, 3]], pce[mine[nm2]["P"]],
    pce[mine[nm2]["P"]] - nbA["SummaryData"][[i, 3]]}], {pair, {{2, "2T", "2T"}, {3, "VMM n=1.88", "VMM"}, {4, "3T", "3T"}, {5, "4T", "4T"}}},
  {i, {pair[[1]]}}, {nm, {pair[[2]]}}, {nm2, {pair[[3]]}}];
nbB = MakeTandemJVPlots["SeriesOptics" -> "SameAsOtherTandems", "n" -> 1.88, "Points" -> 1200,
   "Rz3T" -> 2., "Radd2T" -> 1., "RaddVMMTop" -> 0.5, "RaddVMMBottomString" -> 0.5, "RaddVMMBus" -> 1.,
   "Radd4TTop" -> 1., "Radd4TBottom" -> 1.];
mineB = <|"2T" -> config2T[peR, siR, 1.], "VMM" -> configVMM[peR, siR, 1.88, 0.5, 0.5, 1.],
   "3T" -> config3TL1[peR, siR, 0., 0., 2.], "4T" -> config4T[peR, siR, 1., 1.]|>;
Do[AppendTo[verRows, {"Resistive PCE " <> pair[[2]], nbB["SummaryData"][[pair[[1]], 3]], pce[mineB[pair[[3]]]["P"]],
    pce[mineB[pair[[3]]]["P"]] - nbB["SummaryData"][[pair[[1]], 3]]}],
  {pair, {{2, "2T Radd=1", "2T"}, {3, "VMM top/str 0.5, bus 1", "VMM"}, {4, "3T Rz=2", "3T"}, {5, "4T 1+1", "4T"}}}];
(* thin-optics notebook value for the record *)
nbThin = MakeTandemJVPlots["SeriesOptics" -> "NotebookThin", "n" -> 1.88, "Points" -> 1200];
AppendTo[verRows, {"Notebook 2T with thin-top optics (different generation)", nbThin["SummaryData"][[2, 3]], "", ""}];
AppendTo[verRows, {"Thin-top Jg top / bottom (mA/cm2)", 0.1 jgPerovskitetandemthin, 0.1 jgSitandemthin, ""}];
AppendTo[verRows, {"Max of thin-top EQE (1.1 x EQE)", Max[1.1 RecordPerovskiteEQE[[All, 2]]], "", ""}];
(* recomputed JRad0 at Tref equals the notebook value *)
AppendTo[verRows, {"Generation top (A/m2)", jgPerovskitetandem, peR["g"], peR["g"] - jgPerovskitetandem}];
AppendTo[verRows, {"Generation bottom (A/m2)", jgSitandem, siR["g"], siR["g"] - jgSitandem}];
(* identities *)
AppendTo[verRows, {"3T(Rz=0) - 4T (identity, pp)", 0., "", pce[mine["3T"]["P"]] - pce[mine["4T"]["P"]]}];
AppendTo[verRows, {"3T-L4(m=n) - VMM(n) at R=0 (equivalence, pp)", 0., "",
   pce[config3TL4[peR, siR, 0., 0., 0., nOpt]["P"]] - pce[configVMM[peR, siR, nOpt, 0., 0., 0.]["P"]]}];
(* energy balance at a resistive 3T point: internal power - external = sum R J^2 *)
With[{p = config3TL1[peR, siR, 0.5, 0.5, 2.]},
  AppendTo[verRows, {"3T energy balance residual (W/m2)", 0., "",
    (p["vp"] jf[peR, p["vp"]] + p["vs"] jf[siR, p["vs"]] - p["P"]) - p["Loss_R"]}];
  AppendTo[verRows, {"3T KCL: J_T - J_R - J_Z (A/m2)", 0., "", p["JTop"] - p["JBot"] + p["JZ"]}]];
(* grid convergence *)
{peF, siF} = curvesAt[0., Tref, 1., 1., 0.0005];
resF = evalAll[peF, siF, stNom];
Do[AppendTo[verRows, {"Grid 0.5 mV vs 1 mV (nominal R) " <> c, pce[refNom[c]["P"]], pce[resF[c]["P"]], pce[resF[c]["P"]] - pce[refNom[c]["P"]]}],
  {c, {"2T", "VMM", "3T-L1", "3T-L2", "3T-L4", "4T"}}];
csv["verification", verRows];
log["Verification: ", verRows];

(* ===================================== 3. RESISTANCE SWEEPS (reference) *)
Rgrid = {0., 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1., 1.5, 2., 3., 4., 5.};
resRows = {Join[{"Case", "R_ohm_cm2"}, ("PCE_" <> # & /@ {"2T", "VMM", "3T-L1", "3T-L4", "4T"}),
    ("LossR_W_m2_" <> # & /@ {"2T", "VMM", "3T-L1", "3T-L4", "4T"}), {"JZ_3T_mAcm2", "Jout_VMM_mAcm2", "J_2T_mAcm2"}]};
Do[Module[{st, r},
   st = If[case == "A", <|"Rx" -> R, "Rl" -> 0., "n" -> nOpt, "m" -> m3|>, <|"Rx" -> 0., "Rl" -> R, "n" -> nOpt, "m" -> m3|>];
   r = evalAll[peR, siR, st];
   AppendTo[resRows, Join[{case, R}, pce[r[#]["P"]] & /@ {"2T", "VMM", "3T-L1", "3T-L4", "4T"},
      r[#]["Loss_R"] & /@ {"2T", "VMM", "3T-L1", "3T-L4", "4T"}, {0.1 r["3T-L1"]["JZ"], 0.1 r["VMM"]["Jout"], 0.1 r["2T"]["J"]}]]],
  {case, {"A", "B"}}, {R, Rgrid}];
csv["resistance_sweeps_reference", resRows];
log["Resistance sweeps done"];

(* VMM branch-resistance asymmetry: re-optimised n vs string resistance *)
vmmAsym = {{"Rstr_ohm_cm2", "Rtop_ohm_cm2", "n_opt", "PCE_at_n_opt", "PCE_at_design_n"}};
Do[Module[{nb, pb},
   {nb, pb} = max1D[Function[nn, configVMM[peR, siR, nn, rt, rs, 0.]["P"]], 1.6, 2.2, 60];
   AppendTo[vmmAsym, {rs, rt, nb, pce[pb], pce[configVMM[peR, siR, nOpt, rt, rs, 0.]["P"]]}]],
  {pair, {{0., 0.}, {0.5, 0.}, {1., 0.}, {2., 0.}, {4., 0.}, {0., 0.5}, {0., 1.}, {0., 2.}, {0., 4.}}}, {rs, {pair[[1]]}}, {rt, {pair[[2]]}}];
csv["vmm_branch_resistance", vmmAsym];

(* 3T shared-contact map: Rz x top-generation factor (pure Rz, leads zero) *)
fGrid = Range[0.85, 1.10, 0.025];
RzGrid = {0., 0.5, 1., 2., 3., 5., 7.5, 10., 15., 20.};
mapRows = {{"ftop", "Rz_ohm_cm2", "DeltaJ_mAcm2", "PCE_3T_L1", "PCE_4T_lossless", "Loss_pp", "Rz_DeltaJ2_pp"}};
peByF = Association[Table[f -> (setT[Tref]; buildCurve["Pe", f peR["g"]]), {f, fGrid}]];
Do[Module[{pe = peByF[f], p3, p4},
   p4 = config4T[pe, siR, 0., 0.];
   Do[p3 = config3TL1[pe, siR, 0., 0., rz];
    AppendTo[mapRows, {f, rz, 0.1 p3["JZ"], pce[p3["P"]], pce[p4["P"]], pce[p4["P"]] - pce[p3["P"]],
      100 c2 rz (p4["JTop"] - p4["JBot"])^2/Pin}], {rz, RzGrid}]], {f, fGrid}];
csv["3T_Rz_vs_mismatch_map", mapRows];
log["3T Rz map done"];

(* VMM ratio x mismatch map and 2T mismatch line *)
nMap = Range[1.60, 2.20, 0.02];
vmmMap = {{"ftop", "n", "PCE_VMM", "PCE_4T_lossless", "Loss_pp"}};
Do[Module[{pe = peByF[f], p4},
   p4 = pce[config4T[pe, siR, 0., 0.]["P"]];
   Do[With[{pv = pce[configVMM[pe, siR, nn, 0., 0., 0.]["P"]]}, AppendTo[vmmMap, {f, nn, pv, p4, p4 - pv}]], {nn, nMap}]],
  {f, fGrid}];
csv["VMM_n_vs_mismatch_map", vmmMap];
log["VMM n map done"];

(* =================================== 4. ONE-AT-A-TIME SENSITIVITY SWEEPS *)
sens = {hdr[{"alpha", "APE_eV", "T_C", "ftop", "kNR", "Rstate"}]};
detail = {{"Sweep", "Value", "Rstate", "Config", "PCE", "Loss_vs_4T0_pp", "LossR_W_m2", "JTop_mAcm2", "JBot_mAcm2", "DeltaJ_mAcm2", "VTop_or_Vout", "VBot", "Voc_top", "Voc_bot"}};
runCond[tag_, alpha_, Tc_, f_, kn_] := Module[{pe, si, rI, rN},
   {pe, si} = curvesAt[alpha, Tc + 273.15, f, kn];
   rI = evalAll[pe, si, stIdeal]; rN = evalAll[pe, si, stNom];
   Do[With[{r = rr[[1]], lab = rr[[2]]},
     AppendTo[sens, rowOf[tag, {alpha, ape[alpha], Tc, f, kn, lab}, pe["g"], si["g"], r]];
     Do[If[KeyExistsQ[r, c],
       AppendTo[detail, {tag, Switch[tag, "spectrum", alpha, "temperature", Tc, "mismatch", f, "voltage", kn, _, 0], lab, c,
          pce[r[c]["P"]], pce[r["4T0"]["P"]] - pce[r[c]["P"]], r[c]["Loss_R"],
          0.1 Lookup[r[c], "JTop", Lookup[r[c], "J", 0.]], 0.1 Lookup[r[c], "JBot", Lookup[r[c], "JStr", 0.]],
          0.1 (pe["Jsc"] - si["Jsc"]), Lookup[r[c], "VTZ", Lookup[r[c], "Vout", Lookup[r[c], "VTop", 0.]]],
          Lookup[r[c], "VZR", Lookup[r[c], "VBot", 0.]], pe["Voc"], si["Voc"]}]], {c, Append[cfgOrder, "4T0"]}]],
    {rr, {{rI, "R0"}, {rN, "Rnom"}}}];
   {rI, rN}];
Do[runCond["spectrum", a, 25., 1., 1.], {a, Range[-1., 1., 0.25]}]; log["spectral sweep done"];
Do[runCond["temperature", 0., tc, 1., 1.], {tc, Range[15., 75., 10.]}]; log["temperature sweep done"];
Do[runCond["mismatch", 0., 25., f, 1.], {f, fGrid}]; log["mismatch sweep done"];
Do[runCond["voltage", 0., 25., 1., kn], {kn, {1., 1.*^4, 1.*^5, 3.*^5, 1.*^6, 3.*^6, 1.*^7, 3.*^7}}]; log["voltage sweep done"];
csv["sensitivity_one_at_a_time", sens];
csv["sensitivity_detail", detail];
setT[Tref];

(* VMM n x temperature window *)
vmmT = {{"T_C", "n", "PCE_VMM", "PCE_4T_lossless", "Loss_pp"}};
Do[Module[{pe, si, p4},
   {pe, si} = curvesAt[0., tc + 273.15];
   p4 = pce[config4T[pe, si, 0., 0.]["P"]];
   Do[With[{pv = pce[configVMM[pe, si, nn, 0., 0., 0.]["P"]]}, AppendTo[vmmT, {tc, nn, pv, p4, p4 - pv}]], {nn, nMap}]],
  {tc, Range[15., 75., 10.]}];
csv["VMM_n_vs_temperature_map", vmmT];
log["VMM n-T map done"];

(* 2T window: spectral tilt x top factor *)
twoT = {{"alpha", "APE_eV", "ftop", "PCE_2T", "PCE_4T_lossless", "Loss_pp", "PCE_VMM", "PCE_3T_L1"}};
Do[Module[{gp, gs, si, pe, p4},
   setT[Tref]; {gp, gs} = gens[a, Tref, 1.];
   si = buildCurve["Si", gs];
   Do[pe = buildCurve["Pe", f gp]; p4 = pce[config4T[pe, si, 0., 0.]["P"]];
    AppendTo[twoT, {a, ape[a], f, pce[config2T[pe, si, 0.]["P"]], p4, p4 - pce[config2T[pe, si, 0.]["P"]],
      pce[configVMM[pe, si, nOpt, 0., 0., 0.]["P"]], pce[config3TL1[pe, si, 0., 0., 0.]["P"]]}], {f, Range[0.85, 1.15, 0.025]}]],
  {a, Range[-1., 1., 0.25]}];
csv["2T_alpha_vs_mismatch_map", twoT];
log["2T map done"];

(* ======================================= 5. FACTORIAL ROBUSTNESS SET *)
fac = {hdr[{"alpha", "APE_eV", "T_C", "ftop", "kNR", "Rstate"}]};
Do[Module[{pe, si, rI, rN},
   {pe, si} = curvesAt[a, tc + 273.15, f, 1.];
   rI = evalAll[pe, si, stIdeal]; rN = evalAll[pe, si, stNom];
   AppendTo[fac, rowOf["factorial", {a, ape[a], tc, f, 1., "R0"}, pe["g"], si["g"], rI]];
   AppendTo[fac, rowOf["factorial", {a, ape[a], tc, f, 1., "Rnom"}, pe["g"], si["g"], rN]]],
  {a, {-0.5, 0., 0.5}}, {tc, {25., 45., 65.}}, {f, {0.95, 1., 1.05}}];
csv["factorial_robustness", fac];
setT[Tref];
csv["spectral_tilt_APE", Prepend[Table[{a, ape[a]}, {a, Range[-1., 1., 0.25]}], {"alpha", "APE_350_1050nm_eV"}]];
log["ALL DONE"];
