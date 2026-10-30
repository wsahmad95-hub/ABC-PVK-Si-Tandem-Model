(* Extracts the Input cells of ABC-PVK-Si_Tandem_Model.nb to cells.wl (needs the Mathematica front end).
   Run: wolframscript -file extract_cells.wls *)
work = DirectoryName[$InputFileName];
nb = Import[FileNameJoin[{ParentDirectory[work], "ABC-PVK-Si_Tandem_Model.nb"}], "Notebook"];
cells = Cases[nb, Cell[c_, "Input"|"Code", ___] :> c, Infinity];
conv[b_] := Module[{r = Quiet@UsingFrontEnd[First[FrontEndExecute[FrontEnd`ExportPacket[Cell[b, "Input"], "InputText"]]]]},
  If[StringQ[r], r, ToString[ToExpression[b, StandardForm, HoldComplete], InputForm]]];
strs = MapIndexed[("(* ===== CELL " <> ToString[#2[[1]]] <> " ===== *)\n" <> conv[#1]) &, cells];
note = "(* Input cells of ABC-PVK-Si_Tandem_Model.nb, a modified version of RealisticTandem.nb\n   (Copyright 2017 M. H. Futscher and B. Ehrler; modified 2026 by W. Ahmad). GPL-3.0-or-later, see LICENSE. *)";
strs[[1]] = StringReplace[strs[[1]], " ===== *)\n" -> " ===== *)\n" <> note <> "\n", 1];
Export[FileNameJoin[{work, "cells.wl"}], StringRiffle[strs, "\n\n"], "Text"];
