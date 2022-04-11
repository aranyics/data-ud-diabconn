# data-ud-diabconn
Data for manuscript entitled "Topological Dissimilarities of Hierarchical Resting Networks in Type 2 Diabetes Mellitus and Obesity"


[![DOI](https://zenodo.org/badge/370732223.svg)](https://zenodo.org/badge/latestdoi/370732223)


## Analysis workflow

### Prerequisites

First, fMRI data needs to be preprocessed

### 1_coordAdjust

Independent component analysis-based group-level adjustment of selected region coordinates.

### 2_DCM

Extraction of regional time-series and estimation of effective connectivity parameters with dynamic causal modelling (DCM).

### 3_PEB

Group-level modelling of DCM connectivity with parametric empirical Bayes (PEB)

### 4_graphTheory

Graph theoretical analysis of group-level networks.
