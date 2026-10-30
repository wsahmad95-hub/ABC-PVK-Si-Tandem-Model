# ABC-PVK-Si Tandem Model

Wolfram Mathematica code for the paper

> W. Ahmad and P. Gao, *Interconnection losses and robust design of all-back-contact perovskite/silicon tandem solar cells*.

The code models an all-back-contact (ABC) perovskite/silicon tandem: a quasi-interdigitated back-contact (QIBC) perovskite top cell (1.62 eV, CsFAMA) on a hybrid interdigitated back-contact (HIBC) silicon bottom cell under AM1.5G. It then compares four ways of wiring the subcells:
- two-terminal series (2T);
- voltage-matched module (VMM);
- three-terminal (3T, loading strategies L1–L4);
- four-terminal (4T).

## Origin of the code

This code extends the open-source perovskite/silicon tandem model of M. H. Futscher and B. Ehrler, [Tandem-Solar-Cell-Model](https://github.com/HybridSolarCells/Tandem-Solar-Cell-Model) (`RealisticTandem.nb`, copyright 2017, GPL-3.0). That model is described in:

> M. H. Futscher and B. Ehrler, Modeling the performance limitations and prospects of perovskite/Si tandem solar cells under realistic operating conditions, *ACS Energy Lett.*, 2017, **2**, 2089–2095, [DOI: 10.1021/acsenergylett.7b00596](https://doi.org/10.1021/acsenergylett.7b00596).

The detailed-balance subcell models, the optical coupling and the 2T, VMM and 4T circuits come from that work. The changes made here (2026, W. Ahmad) are:
- the silicon and perovskite subcells were replaced by the HIBC silicon cell and the TCAD-optimized QIBC perovskite cell listed under *Input data*;
- a three-terminal (3T) tandem model and the export of J–V curves and performance tables were added to the notebook;
- the scripts in `analysis/` were added for the wiring analysis of the paper.

The original copyright and license notices are kept at the top of the notebook. If you use this code, please also cite the paper of Futscher and Ehrler.

## Contents

| Path | What it does |
|---|---|
| `ABC-PVK-Si_Tandem_Model.nb` | Main notebook. It builds the detailed-balance models of both subcells from their J–V and EQE data, then calculates the optical coupling and the J–V curves and efficiencies of the 2T, VMM, 3T and 4T tandems. |
| `Input/` | Data read by the notebook (see *Input data* below). |
| `analysis/extract_cells.wls` | Copies the Input cells of the notebook to `analysis/cells.wl`. You only need it after editing the notebook. |
| `analysis/run_model.wls` | Runs `cells.wl` without the front end and saves the model state to `model_state.mx`. |
| `analysis/wiring_analysis.wls` | Wiring analysis of the paper: connection and constraint losses, resistance sweeps, the shared-contact resistance of 3T, voltage-matching ratio sweeps, sensitivity to the operating conditions, the 27 combined conditions and the design windows. Writes `analysis/Results/`. |
| `analysis/summarize.py` | Robustness metrics and design-window limits (`Results/chapter_summary.json`). |
| `analysis/Results/` | The results used in the paper, as CSV and JSON files. |

## Requirements

- Wolfram Mathematica (tested with 14.3; `wolframscript` for the headless scripts)
- Python 3 with NumPy, for `summarize.py` only

## How to run

**Interactively:** open `ABC-PVK-Si_Tandem_Model.nb` in Mathematica and evaluate all cells (*Evaluation > Evaluate Notebook*). The `Input/` folder has to stay next to the notebook.

**The full analysis of the paper**, from the `analysis/` folder:

```bash
wolframscript -file run_model.wls        # subcell and tandem models, about 1.5 min
wolframscript -file wiring_analysis.wls  # wiring analysis, 4-7 min
python summarize.py                      # robustness summary
```

`wiring_analysis.wls` also checks its circuit calculations against the notebook's own tandem routines and writes the comparison to `Results/verification.csv`.

## Input data

| File | Content | Source |
|---|---|---|
| `astmg173.xls` | AM1.5G reference spectrum | ASTM G173-03(2012) |
| `IV-HTIBC-Si.csv`, `EQE-HTIBC-Si.csv` | J–V and EQE of the HIBC silicon cell (27.81% PCE) | G. Wang et al., *Nature*, 2025, **647**, 369–374 |
| `JV-Ideal-1.62.csv`, `EQE-ideal-1.62.csv` | J–V and EQE of the TCAD-optimized QIBC perovskite top cell | W. Ahmad and P. Gao, *Chem. Eng. J.*, 2025, **526**, 171437; the device follows X. Lin et al., *Adv. Energy Mater.*, 2021, **11**, 2100053 |

## Citation

If you use this code, please cite the paper above and the paper of Futscher and Ehrler. The reference to our paper will be completed here once it is published.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. As a modified version of a GPL-3.0 program, it is distributed under the same license. See [LICENSE](LICENSE).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Data taken from the publications listed under *Input data* remain subject to the terms of their sources.
