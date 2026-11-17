# Multi-Omics-Guided Dynamical Systems Modeling of Kupffer Cell–Hepatic Stellate Cell Crosstalk in MASH Fibrosis

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R Version](https://img.shields.io/badge/R-4.6.1-276DC3.svg)](https://www.r-project.org/)

## Overview
This repository contains the signed directed regulatory network topologies, reproducible R analysis pipelines, and processed dynamical modeling outputs accompanying the manuscript:

> **Multi-Omics-Guided Dynamical Systems Modeling Reveals Hysteretic Locking and Dual-Node Therapeutic Vulnerability in Kupffer Cell–Hepatic Stellate Cell Crosstalk During MASH Fibrosis**  
> Kanishtha Kashiv  
> *Department of Bioinformatics and Computational Biology, Parul University, Vadodara, Gujarat, India*

## Repository Structure
- **`network/`**: Signed directed 13-module (15-node), 31-edge Kupffer/SAMac–HSC regulatory circuit in `.topo` (`sRACIPE`) and `.boolnet` (`BoolNet`) formats.
- **`scripts/`**: End-to-end R scripts (`01` through `09`) covering multi-cohort GEO harmonization, single-cell validation, Boolean attractor enumeration, continuous parameter-ensemble ODE modeling (`sRACIPE`, 37,500 trajectories), hysteresis bifurcation sweeps, clinical attractor scoring (n = 416 biopsies), in-silico single/dual-node perturbation screening, and publication figure generation.
- **`results/tables/`**: Harmonized differential expression tables, steady-state attractor centroids, bifurcation trajectories, clinical patient attractor scores, and single/dual-node synergy screens.

## Human Transcriptomic Data Sources (NCBI GEO)
1. **`GSE135251`** — Bulk liver RNA-seq across the MASLD/MASH fibrosis spectrum (n = 216; Govaere et al., *Sci. Transl. Med.* 2020)
2. **`GSE126848`** — Bulk liver RNA-seq across healthy, obese, NAFL, and MASH stages (n = 57; Suppli et al., *Am. J. Physiol.* 2019)
3. **`GSE162694`** — Bulk liver RNA-seq across MASLD/MASH fibrosis stages (n = 143; Pantano et al., *Sci. Rep.* 2021)
4. **`GSE136103`** — Single-cell RNA-seq of healthy vs. cirrhotic human liver non-parenchymal cells (Ramachandran et al., *Nature* 2019)
