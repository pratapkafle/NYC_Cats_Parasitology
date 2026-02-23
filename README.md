# NYC Cats Parasitology

## Reproducible Analysis of Zoonotic Endoparasites in Free-Roaming Cats (NYC)

This repository contains all data and R code required to reproduce the statistical analyses for the study:  
**Zoonotic Endoparasites in Free-Roaming Cats in New York City**

The repository is fully version controlled, uses locked package versions (renv), and runs from a single master script.

## Project Overview

This project performs:  
- Parasite prevalence estimation with 95% confidence intervals  
- Risk factor analysis using Fisher’s Exact Tests  
- Infection intensity modeling using Negative Binomial GLM  
- Generation of publication-ready tables  
- Automatic creation of output directories  

All analyses are conducted in R.

## Repository Structure

```
NYC_Cats_Parasitology/
│
├── run_all.R # Master script (one-command execution)
├── README.md
├── renv.lock # Locked package versions
├── .Rprofile
│
├── Data/
│ └── clean/
│     └── dataset1.csv # Primary dataset
│
├── R/
│ └── FinalDataAnalysisNew.R # Main analysis script
│
├── output/
│ ├── tables/
│ └── figures/
│
└── renv/ # Reproducible package environment
```

## Reproducibility

This project uses:  
- Relative paths via the `here` package  
- Locked package versions via `renv`  
- Automated output folder creation  
- Single-command execution  

This ensures the project runs identically on:  
- A collaborator’s computer  
- A new workstation  
- A server  
- Future systems  

## How to Reproduce the Analysis

### Step 1 — Clone the Repository

```
git clone https://github.com/pratapkafle/NYC_Cats_Parasitology.git
cd NYC_Cats_Parasitology
```

### Step 2 — Restore Exact Package Versions

Open R in the project directory and run:  
```r
renv::restore()
```

This installs the exact package versions used in the study.

### Step 3 — Run Entire Analysis

```r
source("run_all.R")
```

All outputs will be automatically generated in:  
- `output/tables/`  
- `output/figures/`  

## Data

- **File**: `Data/clean/dataset1.csv`  
- Missing values encoded as:  
  - `""`  
  - `"NA"`  
  - `" "` (single space)  
- Infection status defined based on egg per gram (EPG) thresholds as described in the manuscript.

## Statistical Methods

- Exact binomial confidence intervals for prevalence  
- Fisher’s Exact Test for categorical risk factors  
- Negative Binomial GLM for overdispersed egg counts  
- Mann–Whitney U tests when appropriate  
- Automated table export (CSV + formatted output)  

## Output Files

The script automatically generates:  
- `Table1.csv`  
- `Table2.csv`  
- Formatted HTML tables (if enabled)  
- Any figures generated in the analysis  

All results are written to the `output/` directory.

## Computational Environment

- R (version used at time of analysis recorded in `renv.lock`)  
- Packages version-locked via `renv`  
- Session reproducibility ensured through snapshotting  

To inspect package versions:  
```r
renv::status()
```

## Citation

If you use this repository, please cite:  
Kafle P.  
Reproducible Analysis of Zoonotic Endoparasites in Free-Roaming Cats (NYC).  
GitHub Repository.  
(You may optionally archive with Zenodo for DOI minting.)

## Author

Pratap Kafle  
Veterinary Parasitology & One Health Research  

## Reproducibility Statement (For Manuscript Use)

All statistical analyses were conducted in R.  
Code and data required to reproduce all analyses and tables are publicly available at:  
https://github.com/pratapkafle/NYC_Cats_Parasitology  

The computational environment is version-locked using `renv` to ensure reproducibility.
