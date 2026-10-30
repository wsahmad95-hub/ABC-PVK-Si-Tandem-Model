(* Headless run of ABC-PVK-Si_Tandem_Model.nb (its Input cells, extracted to cells.wl): evaluates every cell in order
   and saves the kernel state to model_state.mx for wiring_analysis.wls.  Run: wolframscript -file run_model.wls *)
work = StringReplace[DirectoryName[$InputFileName], "\\" -> "/"];   (* forward slashes: the path is pasted into code *)
base = ParentDirectory[work] <> "/";
work = work <> "/";
Unprotect[NotebookDirectory]; NotebookDirectory[] := base; Protect[NotebookDirectory];
src = Import[work <> "cells.wl", "Text"];

parts = StringDelete[#, StartOfString ~~ DigitCharacter .. ~~ " ===== *)"] & /@
   Select[StringSplit[src, "(* ===== CELL "], StringTrim[#] =!= "" &];
Print["cells: ", Length[parts]];
plotCells = {3, 6, 9, 11, 13, 15, 17, 28, 30, 32, 34, 36, 39, 41, 49, 54, 61, 68};
Do[
  If[MemberQ[plotCells, i], Continue[]];
  code = parts[[i]];
  If[i == 72, code = StringReplace[code, "tandemOutputRoot = Automatic;" -> "tandemOutputRoot = \"" <> work <> "Notebook_Base_Run\";"]];
  t = First@AbsoluteTiming[r = ToExpression[code]];
  If[i == 26, Print["Si single junction: PCE ", SQSi[SiBandgap[T]], " %  Jsc ", jscSi, "  Voc ", vocSi]];
  If[i == 48, Print["Perovskite single junction: PCE ", SQPerovskite[PerovskiteBandgap[T]], " %  Jsc ", jscPerovskite, "  Voc ", vocPerovskite]];
  If[i == 71,
    Print["SQSeries (thin optics): ", SQSeries[SiBandgap[T], PerovskiteBandgap[T]]];
    Print["SQModule n=1: ", SQModule[SiBandgap[T], PerovskiteBandgap[T], 1.]];
    Print["SQThree Rz=0: ", SQThree[SiBandgap[T], PerovskiteBandgap[T], 0.]];
    Print["SQFour: ", SQFour[SiBandgap[T], PerovskiteBandgap[T]]]];
  If[t > 2, Print["cell ", i, " took ", Round[t, 0.1], " s"]],
  {i, Length[parts]}];
Print["PowerofSun ", PowerofSun, " jgSi ", jgSi, " jgPe ", jgPerovskite];
Print["jgPeTandem ", jgPerovskitetandem, " jgSiTandem ", jgSitandem, " thin: ", jgPerovskitetandemthin, " ", jgSitandemthin, " fit ", fit];
Print[tandemRunStatus["JVExport"]["Directory"]];
DumpSave[work <> "model_state.mx", "Global`"];
