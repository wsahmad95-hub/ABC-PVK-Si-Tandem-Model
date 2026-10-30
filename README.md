# ABC-PVK-Si Tandem Model

Wolfram Mathematica code for the paper

> W. Ahmad and P. Gao, *Interconnection losses and robust design of all-back-contact perovskite/silicon tandem solar cells*, submitted to EES Solar.

The code models an all-back-contact (ABC) perovskite/silicon tandem: a quasi-interdigitated back-contact (QIBC) perovskite top cell (1.62 eV, CsFAMA) on a heterojunction interdigitated back-contact (HIBC) silicon bottom cell under AM1.5G. It then compares four ways of wiring the subcells:
- two-terminal series (2T);
- voltage-matched module (VMM);
- three-terminal (3T, loading strategies L1–L4);
- four-terminal (4T).

## Contents

| Path | What it does |
|---|---|
| `ABC-PVK-Si_Tandem_Model.nb` | Main notebook. It builds the detailed-balance models of both subcells from their J–V and EQE data, then calculates the optical coupling and the J–V curves and efficiencies of the 2T, VMM, 3T and 4T tandems. |
| `Input/` | Data read by the notebook (see *Input data* below). |
| `analysis/extract_cells.wls` | Copies the Input cells of the notebook to `analysis/cells.wl`. You only need it after editing the notebook. |
| `analysis/run_model.wls` | Runs `cells.wl` without the front end and saves the model state to `model_state.mx`. |
| `analysis/wiring_analysis.wls` | Wiring analysis of the paper: connection and constraint losses, resistance sweeps, the shared-contact resistance of 3T, voltage-matching ratio sweeps, sensitivity to the operating conditions, the 27 combined conditions and the design windows. Writes `analysis/Results/`. |
| `analysis/make_figures.py` | Draws Figs. 3–11 of the paper from `analysis/Results/` into `analysis/Figures/`. |
| `analysis/summarize.py` | Robustness metrics and design-window limits (`Results/chapter_summary.json`). |
| `analysis/Results/` | The results used in the paper, as CSV and JSON files. |

Figs. 1 and 2 of the paper are schematics and are not produced by the code.

## Requirements

- Wolfram Mathematica (tested with 14.3; `wolframscript` for the headless scripts)
- Python 3 with NumPy and Matplotlib, for the figures only (tested with Python 3.13, NumPy 2.3, Matplotlib 3.10). The figures use the Times New Roman font.

## How to run

**Interactively:** open `ABC-PVK-Si_Tandem_Model.nb` in Mathematica and evaluate all cells (*Evaluation > Evaluate Notebook*). The `Input/` folder has to stay next to the notebook.

**The full analysis of the paper**, from the `analysis/` folder:

```bash
wolframscript -file run_model.wls        # subcell and tandem models, about 1.5 min
wolframscript -file wiring_analysis.wls  # wiring analysis, 4-7 min
python make_figures.py                   # Figs. 3-11
python summarize.py                      # robustness summary
```

`wiring_analysis.wls` also checks its circuit calculations against the notebook's own tandem routines and writes the comparison to `Results/verification.csv` (Table S3 of the ESI).

## Figures of the paper

| Paper | File written by `make_figures.py` |
|---|---|
| Fig. 3 | `Fig3_reference_JV.png` |
| Fig. 4 | `Fig4_loss_decomposition.png` |
| Fig. 5 | `Fig5_resistance_sweeps.png` |
| Fig. 6 | `Fig6_3T_shared_contact.png` |
| Fig. 7 | `Fig7_ratio_sweeps.png` |
| Fig. 8 | `Fig8_sensitivity_constraint_loss.png` |
| Fig. 9 | `Fig9_sensitivity_nominal_R.png` |
| Fig. 10 | `Fig10_factorial_robustness.png` |
| Fig. 11 | `Fig11_design_windows.png` |

## Input data

| File | Content | Source |
|---|---|---|
| `astmg173.xls` | AM1.5G reference spectrum | ASTM G173-03(2012) |
| `IV-HTIBC-Si.csv`, `EQE-HTIBC-Si.csv` | J–V and EQE of the 27.4% HIBC silicon cell | G. Wang et al., *Nature*, 2025, **647**, 369–374 |
| `JV-Ideal-1.62.csv`, `EQE-ideal-1.62.csv` | J–V and EQE of the TCAD-optimized QIBC perovskite top cell | W. Ahmad and P. Gao, *Chem. Eng. J.*, 2025, **526**, 171437; the device follows X. Lin et al., *Adv. Energy Mater.*, 2021, **11**, 2100053 |

## Citation

If you use this code, please cite the paper above. The reference will be updated here once the paper is published.

## License

MIT, see [LICENSE](LICENSE). Data taken from the publications listed under *Input data* remain subject to the terms of their sources.
