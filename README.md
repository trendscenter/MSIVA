# MSIVA: Multimodal Subspace Independent Vector Analysis

This repository contains a MATLAB implementation for Multimodal Subspace Independent Vector Analysis (MSIVA).

## Simulation

Simulation code can be found in [`simulation.m`](scripts/SIVA/simulation.m).

MSIVA initialization workflow is implemented in [`run_mgpca_ica.m`](scripts/SIVA/run_mgpca_ica.m). Unimodal initialization workflow is implemented in [`run_pca_ica.m`](scripts/SIVA/run_pca_ica.m). Multimodal initialization workflow is implemented in [`run_mgpca_gica.m`](scripts/SIVA/run_mgpca_gica.m).

## Neuroimaging Experiment

Neuroimaging experiment code can be found in [`experiment.m`](scripts/SIVA/experiment.m).

## Visualization

Code for the ISBI 2023 paper figures can be found in [`figures/ISBI2023SIVA`](figures/ISBI2023SIVA).

**Fig. 1**: [`plot_subspace_struct.ipynb`](figures/ISBI2023SIVA/plot_subspace_struct.ipynb)

**Fig. 2**: [`plot_sim.ipynb`](figures/ISBI2023SIVA/plot_sim.ipynb)

**Fig. 3**: [`plot_sim.ipynb`](figures/ISBI2023SIVA/plot_sim.ipynb)

**Fig. 4**: [`plot_img.ipynb`](figures/ISBI2023SIVA/plot_img.ipynb)

**Fig. 5**: [`plot_img.ipynb`](figures/ISBI2023SIVA/plot_img.ipynb)

**Fig. 6**: [`dualmap.m`](figures/ISBI2023SIVA/dualmap.m)

Code for the bioRxiv preprint figures can be found in [`figures/Journal2024MSIVA`](figures/Journal2024MSIVA).

**Fig. 1**: [`plot_subspace_struct.ipynb`](figures/Journal2024MSIVA/plot_subspace_struct.ipynb)

**Fig. 3**: [`plot_sim.ipynb`](figures/Journal2024MSIVA/plot_sim.ipynb)

**Fig. 4**: [`plot_sim.ipynb`](figures/Journal2024MSIVA/plot_sim.ipynb)

**Fig. 5**: [`plot_img_ukb.ipynb`](figures/Journal2024MSIVA/plot_img_ukb.ipynb)

**Fig. 6**: [`plot_img_sz.ipynb`](figures/Journal2024MSIVA/plot_img_sz.ipynb)

**Fig. 7**: [`plot_img_ukb.ipynb`](figures/Journal2024MSIVA/plot_img_ukb.ipynb)

**Fig. 8**: [`plot_img_sz.ipynb`](figures/Journal2024MSIVA/plot_img_sz.ipynb)

**Figs. 9 & 10**: [`dualcodeImage_AY_geomedian.m`](figures/Journal2024MSIVA/dualcodeImage_AY_geomedian.m)

**Fig. 11**: [`plot_sig_voxel.ipynb`](figures/Journal2024MSIVA/plot_sig_voxel.ipynb)

**Fig. 12**: (1) Run [`age_delta.m`](figures/Journal2024MSIVA/age_delta.m) to compute brain-age delta. (2) Run [`compute_geometric_median.py`](figures/Journal2024MSIVA/compute_geometric_median.py) to compute geometric median of brain-age delta. (3) Run [`phenotype_map.py`](figures/Journal2024MSIVA/phenotype_map.py) to compute spatial correlation between brain-age delta and phenotype variable. (4) Use [`dualcodeImage_beta1.m`](figures/Journal2024MSIVA/dualcodeImage_beta1.m), [`dualcodeImage_delta2p_std.m`](figures/Journal2024MSIVA/dualcodeImage_delta2p_std.m), [`dualcodeImage_delta2p_geomedian.m`](figures/Journal2024MSIVA/dualcodeImage_delta2p_geomedian.m), and [`dualcodeImage_phenotype.m`](figures/Journal2024MSIVA/dualcodeImage_phenotype.m) to plot the dual-coded maps.

**Fig. 13**: [`dualcodeImage_beta1.m`](figures/Journal2024MSIVA/dualcodeImage_beta1.m)

**Fig. 14**: [`plot_img_ukb_rdc.ipynb`](figures/Journal2024MSIVA/plot_img_ukb_rdc.ipynb)

**Fig. 15**: [`plot_img_sz_rdc.ipynb`](figures/Journal2024MSIVA/plot_img_sz_rdc.ipynb)

**Figs. 16 & 17**: [`dualcodeImage_AY_geomedian.m`](figures/Journal2024MSIVA/dualcodeImage_AY_geomedian.m)

**Fig. 18**: [`plot_num_crossmodal_voxel.ipynb`](figures/Journal2024MSIVA/plot_num_crossmodal_voxel.ipynb)

**Figs. 19 & 20**: [`compare_mmiva_msiva_ukb.ipynb`](figures/Journal2024MSIVA/compare_mmiva_msiva_ukb.ipynb)

**Figs. 21 & 22**: [`compare_mmiva_msiva_sz.ipynb`](figures/Journal2024MSIVA/compare_mmiva_msiva_sz.ipynb)

## Prerequisites

[Group ICA Of fMRI Toolbox (GIFT)](http://trendscenter.org/software/gift/)

[Multidataset Independent Subspace Analysis (MISA)](https://github.com/rsilva8/MISA)

## References
If you find this repository useful, please consider citing at least one of the following papers:
```
@article {li2024multimodal,
	author = {Li, Xinhui and Kochunov, Peter and Adali, Tulay and Silva, Rogers F and Calhoun, Vince},
	title = {Multimodal subspace independent vector analysis effectively captures the latent relationships between brain structure and function},
	elocation-id = {2023.09.17.558092},
	year = {2024},
	doi = {10.1101/2023.09.17.558092},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/10/22/2023.09.17.558092},
	eprint = {https://www.biorxiv.org/content/early/2024/10/22/2023.09.17.558092.full.pdf},
	journal = {bioRxiv}
}

@inproceedings{li2023multimodal,
  title={Multimodal subspace independent vector analysis better captures hidden relationships in multimodal neuroimaging data},
  author={Li, Xinhui and Adali, Tulay and Silva, Rogers F and Calhoun, Vince D},
  booktitle={2023 IEEE 20th International Symposium on Biomedical Imaging (ISBI)},
  pages={1--5},
  year={2023},
  organization={IEEE}
}

@article{silva2020multidataset,
  title={Multidataset independent subspace analysis with application to multimodal fusion},
  author={Silva, Rogers F and Plis, Sergey M and Adal{\i}, T{\"u}lay and Pattichis, Marios S and Calhoun, Vince D},
  journal={IEEE Transactions on Image Processing},
  volume={30},
  pages={588--602},
  year={2020},
  publisher={IEEE}
}
```