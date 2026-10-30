(* Extracts the Input cells of ABC-PVK-Si_Tandem_Model.nb to cells.wl (needs the Mathematica front end).
   Run: wolframscript -file extract_cells.wls *)
work = DirectoryName[$InputFileName];
nb = Import[FileNameJoin[{ParentDirectory[work], "ABC-PVK-Si_Tandem_Model.nb"}], "Notebook"];
cells = Cases[nb, Cell[c_, "Input"|"Code", ___] :> c, Infinity];
conv[b_] := Module[{r = Quiet@UsingFrontEnd[First[FrontEndExecute[FrontEnd`ExportPacket[Cell[b, "Input"], "InputText"]]]]},
  If[StringQ[r], r, ToString[ToExpression[b, StandardForm, HoldComplete], InputForm]]];
strs = MapIndexed[("(* ===== CELL " <> ToString[#2[[1]]] <> " ===== *)\n" <> conv[#1]) &, cells];
Export[FileNameJoin[{work, "cells.wl"}], StringRiffle[strs, "\n\n"], "Text"];
