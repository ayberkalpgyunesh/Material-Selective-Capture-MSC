# R scripts

Statistics and plotting for the manuscript figures. Each script takes the extracted `*_r_input/` folder from the matching notebook bundle as its input directory (command-line argument or `INPUT_DIR` at the top of the script). No local paths are hard-coded.

```
Rscript figure2_final.r /path/to/05_r_input
```

| Script | Manuscript figures | Input folder (from notebook) |
|---|---|---|
| `figure2_final.r` | Fig. 2, S1 Fig | `05_r_input/` from `01_Fig2_Fig3_S1_S5_docking_and_variants.ipynb` |
| `figure3_final.r` | Fig. 3, S5 Fig | `05_r_input/` from `01_Fig2_Fig3_S1_S5_docking_and_variants.ipynb` |
| `figure4_final.r` | Fig. 4 (plotting; Wilcoxon statistics are computed in the notebook) | `04_r_input/` from `02_Fig4_hysteresis.ipynb` |
| `figure5_final.r` | Fig. 5 (seed-clustered logistic models, Holm, TOST equivalence) | `04_r_input/` from `03_Fig5_anti_tautology.ipynb` |
| `MSC_Figure6_and_FigureS4_panels.r` | Fig. 6, S4 Fig | `10_r_input/` folders from `04a`, `04b`, `04c` |
| `MSC_SupplementaryFigureS_material_dynamics.r` | S2 Fig | `10_r_input/` from `05_S2_material_dynamics.ipynb` |

S3 Fig is produced directly in Python by `06_S3_sensitivity.ipynb`.

Note: some script headers use working names from earlier drafts ("Supplementary Figure S9" = S4 Fig, "Supplementary Figure S1" in the material-dynamics script = S2 Fig). The table above gives the manuscript numbering.

Main packages: ggplot2, patchwork, dplyr, data.table, arrow, fixest.
