# TNBC Spatial LDA

This repository contains R and Stan code for exploring topic-model-based representations of the triple-negative breast cancer (TNBC) tumor microenvironment using spatially resolved single-cell imaging data.

The project compares conventional Latent Dirichlet Allocation (LDA) with a spatially informed LDA workflow. In this setting, documents are defined either at the image/sample level or at the spatial-tile level, and words correspond to observed cell phenotypes. The spatial LDA workflow introduces spatial information into the topic model by using tile locations and topic centers to construct a spatially informed prior for document-topic proportions.

## Abstract

In cancer research, characterizing the tumor microenvironment (TME) is essential for understanding tumor behavior and patient outcomes. Advanced spatial imaging technologies, such as Multiplexed Ion Beam Imaging by Time of Flight (MIBI-TOF), provide detailed single-cell spatial information about the TME. However, conventional Latent Dirichlet Allocation (LDA) does not explicitly account for spatial proximity and may identify topics that are weakly separated or difficult to interpret in spatially complex tissue data.

This project introduces a spatial Latent Dirichlet Allocation (spatial LDA) workflow for modeling cellular communities in triple-negative breast cancer (TNBC) tissue images. In this approach, spatially localized documents are constructed from image samples and cell coordinates, with each document represented as a mixture of latent topics. Topic structure is informed by spatial proximity among documents, allowing the model to better reflect local organization within the TME. By comparing spatial LDA with conventional LDA, this repository evaluates whether incorporating spatial information improves the interpretability and spatial coherence of inferred tumor microenvironment topics.

## Research Motivation

The tumor microenvironment (TME) is spatially organized: immune, stromal, and tumor cell populations do not occur independently of location. Standard LDA can identify latent cellular communities from cell-type counts, but it does not directly account for spatial proximity. This repository explores whether incorporating spatial information into topic modeling can produce more spatially coherent and interpretable representations of TNBC tissue structure.

The main goals are to:

- preprocess TNBC multiplexed imaging data into spatial and cell-type representations;
- fit conventional LDA models to cell phenotype counts;
- construct spatial tiles from cell coordinates;
- estimate spatially informed topic distributions using Stan;
- compare topic structures, spatial patterns, and model diagnostics between standard LDA and spatial LDA.

## Repository Structure

```text
TNBC_Spatial_LDA/
├── README.md
├── Rstan/
│   ├── lda.stan
│   └── slda.stan
└── R_notebooks/
    ├── 01_TNBC_Exploratory_LC.Rmd
    ├── 02_TNBC_spe_LC.Rmd
    ├── 03_TNBC_LMM.Rmd
    ├── 04_MB_Xeno_LDA/
    ├── 05_LDA_multiChains/
    ├── 06_LDA_scripts/
    │   ├── 06_LDA_analysis.R
    │   ├── 06_lda.stan
    │   ├── 06_alignmentMatrix.R
    │   ├── 06_betaAligned.R
    │   ├── 06_thetaAligned.R
    │   ├── 06_diagnosticsPlot.R
    │   └── model-assessment and plotting utilities
    ├── 07_TNBC_Diagnostics.Rmd
    ├── 08_TNBC_topic_proportion.Rmd
    ├── 09_TNBC_topic_LMM.Rmd
    ├── 10_spatial clustering.Rmd
    ├── 11_tile_dataframe.Rmd
    ├── 12_Spatial_LDA/
    │   ├── 12_spatial_LDA.Rmd
    │   ├── 12_spatial_LDA_1000tiles_Patient 4.Rmd
    │   ├── 12_spatial_LDA_1000tiles_Patient 12.Rmd
    │   └── multi-chain patient-specific spatial LDA notebooks
    ├── 13_slda_diagnostics/
    │   ├── 13_Patient 4_slda_diagnostics.Rmd
    │   ├── 13_Patient 12_slda_diagnostics copy.Rmd
    │   └── 13_spatial_lda_assesment.Rmd
    ├── TNBC_Tesselation_V2.Rmd
    └── TNBC_spiat.Rmd
```

## Methodological Overview

### 1. Data Representation

The analysis starts from spatially resolved TNBC cell-level data. Each cell has:

- a sample or patient identifier;
- a cell phenotype or cell-type label;
- spatial coordinates;
- additional cell-level attributes such as cell size.

The data are converted into a `SpatialExperiment` object for downstream spatial analysis and visualization.

### 2. Standard LDA

The conventional LDA model represents each sample as a mixture of latent topics, where each topic is a distribution over cell phenotypes.

In this model:

- documents are samples or spatial units;
- words are cell phenotypes;
- `theta` represents document-topic proportions;
- `beta` or `phi` represents topic-cell-type distributions.

The main standard LDA workflow is implemented in:

```text
R_notebooks/06_LDA_scripts/06_LDA_analysis.R
R_notebooks/06_LDA_scripts/06_lda.stan
R_notebooks/08_TNBC_topic_proportion.Rmd
```

The standard LDA outputs include:

- posterior samples of topic proportions;
- topic-cell-type composition heatmaps;
- sample-topic distribution heatmaps;
- spatial visualization of topic assignments or topic probabilities;
- posterior diagnostics, including R-hat and effective sample size.

### 3. Spatial Tile Construction

To incorporate spatial organization, individual cells are grouped into spatial tiles. Each tile is treated as a document, and the cell types observed within that tile form the document words.

The tile-level workflow includes:

1. extracting cell coordinates from the spatial object;
2. clustering cells into spatial tiles;
3. computing tile centroids;
4. counting cell phenotypes within each tile;
5. converting tile-level cell-type counts into long-format LDA input.

This workflow is mainly developed in:

```text
R_notebooks/11_tile_dataframe.Rmd
R_notebooks/12_Spatial_LDA/
```

### 4. Spatial LDA

The spatial LDA model extends the standard LDA structure by allowing the document-topic prior to depend on spatial information. In the current implementation, spatial influence is computed from the distance between tile centroids and estimated topic centers.

The spatial influence matrix is then passed to Stan and used as the Dirichlet prior for the document-topic distributions.

The main spatial LDA files are:

```text
Rstan/slda.stan
R_notebooks/12_Spatial_LDA/
R_notebooks/13_slda_diagnostics/
```

The spatial LDA outputs include:

- tile-level topic proportions;
- spatial topic probability maps;
- cell-level visualizations mapped back from tile-level estimates;
- MCMC diagnostics for `theta` and `phi`;
- posterior predictive checks comparing observed and simulated cell-type counts.

## Statistical Models

### Standard LDA

The standard LDA model assumes:

```text
theta_d ~ Dirichlet(alpha)
phi_k   ~ Dirichlet(gamma)
word_n  ~ Categorical(theta_d * phi)
```

where:

- `d` indexes documents;
- `k` indexes topics;
- `theta_d` is the topic mixture for document `d`;
- `phi_k` is the cell-type distribution for topic `k`.

### Spatial LDA

The spatial LDA model keeps the same topic-model structure but replaces the fixed Dirichlet prior for document-topic proportions with a spatially informed prior:

```text
theta_d ~ Dirichlet(spatial_influence_d)
phi_k   ~ Dirichlet(gamma)
word_n  ~ Categorical(theta_d * phi)
```

where `spatial_influence_d` summarizes the proximity of document `d` to estimated topic centers.

## Data Requirements

The raw data and some intermediate `.rds` / `.RData` objects are not included in this repository. To run the full workflow, you need to provide the required TNBC imaging data and processed objects in the expected directory structure.

The notebooks currently refer to files such as:

```text
Data/MIBI-TNBC_count_output.csv
Output/RData/02_TNBC_spe_LC.rds
Output/RData/05_LDA_multiChains/
```

Some notebooks also contain local absolute paths. These should be replaced with project-relative paths using `here::here()` before running the analysis on a new machine.

Recommended project structure:

```text
TNBC_Spatial_LDA/
├── Data/
│   └── MIBI-TNBC_count_output.csv
├── Output/
│   └── RData/
│       ├── 02_TNBC_spe_LC.rds
│       └── 05_LDA_multiChains/
├── R_notebooks/
└── Rstan/
```

## Software Requirements

The analysis is written primarily in R, with Bayesian models implemented in Stan through `rstan`.

Main R packages used across the notebooks include:

```r
rstan
dplyr
tidyr
tibble
purrr
ggplot2
readr
magrittr
stringr
reshape2
abind
here
SpatialExperiment
EBImage
spatstat
sp
spdep
adespatial
randomcoloR
```

Because some spatial R packages have changed or been retired over time, it is recommended to record the working package versions using `renv` or another reproducible environment manager.

## Installation

Clone the repository:

```bash
git clone https://github.com/cuix23/TNBC_Spatial_LDA.git
cd TNBC_Spatial_LDA
```

Install the required R packages:

```r
install.packages(c(
  "rstan", "dplyr", "tidyr", "tibble", "purrr",
  "ggplot2", "readr", "magrittr", "stringr",
  "reshape2", "abind", "here", "sp", "spdep",
  "adespatial", "randomcoloR"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "SpatialExperiment",
  "EBImage"
))
```

You may also need a working C++ toolchain for `rstan`.

## Suggested Workflow

A typical workflow is:

1. **Exploratory analysis**

   ```text
   R_notebooks/01_TNBC_Exploratory_LC.Rmd
   ```

2. **Create or load the SpatialExperiment object**

   ```text
   R_notebooks/02_TNBC_spe_LC.Rmd
   ```

3. **Fit standard LDA and summarize topic proportions**

   ```text
   R_notebooks/06_LDA_scripts/06_LDA_analysis.R
   R_notebooks/08_TNBC_topic_proportion.Rmd
   ```

4. **Construct spatial tiles**

   ```text
   R_notebooks/11_tile_dataframe.Rmd
   ```

5. **Run spatial LDA**

   ```text
   R_notebooks/12_Spatial_LDA/
   Rstan/slda.stan
   ```

6. **Evaluate diagnostics and posterior predictive checks**

   ```text
   R_notebooks/13_slda_diagnostics/
   ```

## Outputs

The workflow produces several types of outputs:

- document-topic posterior summaries;
- topic-cell-type composition summaries;
- topic assignment visualizations;
- spatial maps of topic probabilities;
- tile-level and patient-level topic summaries;
- R-hat and effective sample size diagnostics;
- posterior predictive checking plots.
