# MSIVA: Multimodal Subspace Independent Vector Analysis

This repository contains a MATLAB implementation for Multimodal Subspace Independent Vector Analysis (MSIVA), which is built upon [Multidataset Independent Subspace Analysis (MISA)](https://github.com/rsilva8/MISA). 

## Simulation

Simulation code can be found in [`simulation.m`](scripts/SIVA/simulation.m).

Unimodal analysis is implemented in [`run_pca_ica.m`](scripts/SIVA/run_pca_ica.m). MSIVA is implemented in [`run_mgpca_ica.m`](scripts/SIVA/run_mgpca_ica.m).

## Neuroimaging Experiment

Neuroimaging experiment code can be found in [`experiment.m`](scripts/SIVA/experiment.m).

## Visualization

Visualization code for ISBI 2023 paper figures can be found in [`figures/ISBI2023SIVA`](figures/ISBI2023SIVA).

**Fig. 1**: [`plot_subspace_struct.ipynb`](figures/ISBI2023SIVA/plot_subspace_struct.ipynb)

**Fig. 2**: [`plot_sim.ipynb`](figures/ISBI2023SIVA/plot_sim.ipynb)

**Fig. 3**: [`plot_sim.ipynb`](figures/ISBI2023SIVA/plot_sim.ipynb)

**Fig. 4**: [`plot_img.ipynb`](figures/ISBI2023SIVA/plot_img.ipynb)

**Fig. 5**: [`plot_img.ipynb`](figures/ISBI2023SIVA/plot_img.ipynb)

**Fig. 6**: [`dualmap.m`](figures/ISBI2023SIVA/dualmap.m)

## Dependency

[GIFT](http://trendscenter.org/software/gift/) toolbox is required to run MSIVA.