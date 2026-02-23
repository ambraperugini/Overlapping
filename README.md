# README: How do my distributions differ?  
## Significance testing for the Overlapping Index using Permutation Test

### Project Overview
This repository contains all materials required to reproduce the simulations, analyses, and empirical example presented in the article:

The repository is fully reproducible and includes:
- executable article sources (`.Rnw`) used to generate the manuscript;
- simulation code and pre-generated datasets;
- versions of the manuscript across the review process;
- the real dataset used as an applied example in the paper.

---

## Repository Structure

### 1. RNW PAPER (Executable Manuscripts)
The complete history of executable versions of the manuscript written in Sweave (`.Rnw`):
 `.Rnw` file reproduces the full paper, including:
- statistical analyses,
- simulations,
- figures and tables,
- embedded code used to generate results.

Compiling these files produces the corresponding PDF versions of the manuscript.

### 2. Simulation Code
File: `simulation_code.R`

This script contains the full simulation framework used in the article.  
Main components include:

- generation of skewed datasets using the skew-normal distribution (`sn` package);
- implementation of permutation tests (overlapping index, mean-based tests, variance tests);
- configurable simulation parameters:
  - sample size (`n`)
  - skewness (`alpha`)
  - scale (`sigma`)
  - location differences.

Running this script allows users to recreate or extend the simulation study.

### 3. Simulated Data
File: `R02_sim07.rda`

Pre-generated datasets covering all parameter combinations analyzed in the article.

Purpose:
- allows immediate reproduction of results without rerunning computationally intensive simulations.

Structure:
- each row corresponds to one simulation condition;
- columns contain descriptive statistics and test outcomes (means, variances, skewness, test statistics, p-values).

### 4. Reproducible Article Source
File: `Article_overlapping_Main.Rnw`

This Sweave document reproduces the complete article workflow:

- methodological explanations,
- simulation analyses,
- empirical example,
- automatic generation of figures and tables.

#### Compilation
1. Open the `.Rnw` file in RStudio.
2. Compile using **Knit** or a LaTeX-compatible workflow.
3. The output PDF reproduces the article exactly.

### 5. Empirical Dataset (Applied Example)
File: `EngTurk.csv`

This dataset corresponds to the real-data example discussed in the main draft (see references in the manuscript).  
It is used to demonstrate the applied workflow of the overlapping approach on experimentally collected data.
Oksuz, D. C., & Rebuschat, P. (2024, Mar). Collocational processing in typologically different languages,
english and turkish. OSF. Retrieved from osf.io/muwjz

The dataset allows readers to:
- reproduce the empirical example shown in the paper,
- rerun permutation analyses,
- replicate figures and reported statistics directly from raw data.

## Dependencies

Required R packages:

- `sn` — skew-normal data generation
- `overlapping` — overlap statistics
- `parallel` — parallel computation
- `ggplot2` — visualization
- `knitr` — compilation of `.Rnw` documents

### Installation
```r
install.packages(c("sn", "overlapping", "parallel", "ggplot2", "knitr"))
