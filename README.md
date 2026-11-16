# Multi-omic signatures of traffic-related air pollution in London

Analysis code for:

> Cheng SL, Hedges M, Keski-Rahkonen P, *et al.* (2024). **Multiomic signatures of traffic-related air pollution in London reveal potential short-term perturbations in gut microbiome-related pathways.** *Environmental Science & Technology*, 58(20), 8771–8782. [doi:10.1021/acs.est.3c09148](https://doi.org/10.1021/acs.est.3c09148) ([correction](https://doi.org/10.1021/acs.est.4c09132))

## Overview

The Oxford Street II study is a randomised crossover study. Participants walked for two hours on Oxford Street, a busy London street open only to diesel buses and taxis, and in nearby traffic-free Hyde Park. Blood was sampled 2 h before and 2 h and 24 h after each walk. Personal exposure to PM<sub>10</sub>, PM<sub>2.5</sub>, NO<sub>2</sub>, black carbon (BC) and particle number count (PCNT) was measured during the walks. The analyses include 50 participants (healthy, COPD or ischaemic heart disease) with serum untargeted metabolomics (6,040 features, UHPLC-QTOF-MS) and, in a subset, blood transcriptomics (30,923 transcripts).

| Folder | Scripts | Analysis |
|---|---|---|
| `01-Preprocessing` | 01–05 | QRILC imputation, data cleaning, removal of technical variation with linear mixed models, probe annotation, PCA |
| `02-Univariate` | 06–08 | Multivariate normal models (generalised least squares for repeated measures) for each feature and exposure, adjusted for age, sex, BMI and health group; feature selection and plots |
| `03-Clustering` | 09–10 | Stability-calibrated consensus clustering ([`sharp`](https://cran.r-project.org/package=sharp)) of the selected metabolic features, summarised by their medoids |
| `04-Pathways` | 11 | Mummichog pathway analysis (run in [MetaboAnalyst](https://www.metaboanalyst.ca)) |
| `05-Networks` | 12–13 | Stability-calibrated graphical LASSO networks of metabolic features and transcripts at each time point |
| `06-Sensitivity` | 14–16 | Site, temperature, humidity and noise, and outlying samples |

Shared functions and settings (significance thresholds, number of cores) are in `functions.R`. Each script starts with a header listing its inputs and outputs.

## Running the code

The participant-level data are not included in this repository. The scripts expect this repository to be cloned as `Scripts/` inside a project folder that holds the data:

```
project/
├── Data/                      raw data
├── New_data/                  processed data (created by 01-Preprocessing)
├── Results/                   results and figures (created by the scripts)
├── Supplementary_Materials/   SI_Table_5.xlsx (metabolite annotation, SI Table S5 of the paper)
└── Scripts/                   this repository
```

Run the scripts in numerical order from the project folder, or set `OXFORD_DIR` to its path. `N_CORES` sets the number of cores for the association models (default 5). Two steps are done outside R: the feature annotation (`Results/metab_annotation_1309.csv`) and the Mummichog analysis in MetaboAnalyst (see the header of `04-Pathways/11-Mummichog_pathways.R`).

## Requirements

R (≥ 4.1) with the following packages:

- **CRAN:** `dplyr`, `tidyr`, `stringr`, `tibble`, `readxl`, `ggplot2`, `patchwork`, `ggrepel`, `RColorBrewer`, `ComplexUpset`, `pheatmap`, `ggvenn`, `circlize`, `imputeLCMD`, `lme4`, `nlme`, `FactoMineR`, `factoextra`, `sharp`, `igraph`
- **Bioconductor:** `GEOquery`, `ComplexHeatmap`; `pcaMethods` and `impute` are needed by `imputeLCMD`

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "ComplexHeatmap", "pcaMethods", "impute"))
install.packages(c(
  "dplyr", "tidyr", "stringr", "tibble", "readxl", "ggplot2", "patchwork", "ggrepel",
  "RColorBrewer", "ComplexUpset", "pheatmap", "ggvenn", "circlize", "imputeLCMD", "lme4",
  "nlme", "FactoMineR", "factoextra", "sharp", "igraph"
))
```

## Citation

If you use this code, please cite the paper:

```bibtex
@article{cheng2024trap,
  title   = {Multiomic Signatures of Traffic-Related Air Pollution in London Reveal Potential
             Short-Term Perturbations in Gut Microbiome-Related Pathways},
  author  = {Cheng, Sibo Lucas and Hedges, Michael and Keski-Rahkonen, Pekka
             and Chatziioannou, Anastasia Chrysovalantou and Scalbert, Augustin and Chung, Kian Fan
             and Sinharay, Rudy and Green, David C. and de Kok, Theo M. C. M. and Vlaanderen, Jelle
             and Kyrtopoulos, Soterios A. and Kelly, Frank and Portengen, L{\"u}tzen and Vineis, Paolo
             and Vermeulen, Roel C. H. and Chadeau-Hyam, Marc and Dagnino, Sonia},
  journal = {Environmental Science \& Technology},
  volume  = {58},
  number  = {20},
  pages   = {8771--8782},
  year    = {2024},
  doi     = {10.1021/acs.est.3c09148}
}
```

## Licence

MIT, see [LICENSE](LICENSE).

## Author

Sibo Lucas Cheng ([ORCID 0009-0009-8930-6148](https://orcid.org/0009-0009-8930-6148)), Department of Epidemiology and Biostatistics, School of Public Health, Imperial College London.
