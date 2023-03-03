# Example scripts
These are scripts to run various examples included in the companion MISA publication.

## Synthetic data examples
### `ICA1` (`V > N` case):
This example assesss the effect of additive sensor noise and condition number in a moderately large independenc component analysis (ICA) problem with *rectangular* mixing matrix `A` and fairly low number of observations `N`. The experiment was setup with `M = 1` dataset containing `C = C_1 = 75` sources organized into `K = C = 75` one-dimensional subspaces, and `N = 3500` observations sampled independently from a Laplace distribution. In each of the ten instances of this experiment, a new, unique (`V x C`) rectangular mixing matrix `A` (`V = 8000`) was randomly generated. For each instance, ten runs were performed, each with a different random row-orthogonal `W0` initialization (thus, 100 runs per experimental condition). The experiments were broken into two parts:

- Part a) Gaussian white noise `e` was added to the mixtures to yield a `SNR = 3dB`, while varying the condition number: `cond(A) = {1,3,7,15}`.
- Part b) the condition number was fixed at `cond(A) = 7`, while Gaussian white noise `e` was added to the mixtures with varying `SNR = {30, 10, 3, 0.4, 0.004} dB`. These SNR values correspond to `signal power : noise power ratios` of `[999:1, 99:1, 1:1, 1:99, 1:999]`.

Here, pseudoinverse reconstructino error (PRE) was employed for data reduction, followed by either Infomax or MISA for final separation.

`Optimization parameters`: For the data reduction, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of `W` were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 50000, the maximum number of iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for `L-BFGS-B` was set to 5. `TolX` and `TolFun`
were set as described in Section II-B of the supplemental material for the acompanying publication, with `b = 80` for part a), and `b = 4` for part b). For MISA, we utilized the same parameters, except: `TolFun = 10^-4`, `TolX = 10^-9`, the initial barrier parameter was set to 1, and the number of past gradients utilized for `L-BFGS-B` was set to 20. For Infomax, we utilized all the defaults in [GIFT] and `sphering` was set to `on`.

[GIFT]: http://trendscenter.org/software/gift/

### `IVA1` (`V < N`, `V = C` case):
In this experiment, we want to assess the performance in an independent vector analysis (IVA) problem when no data reduction is required (i.e., `V = C`), noise is absent, and the number of observations `N` is abundant. The experiment was setup with `M = 10` datasets, each containing `C_m = 16` sources for a total of `K = Cm = 16` ten-dimensional subspaces, and `N = 32968` observations sampled independently from a multivariate Laplace distribution.

The selected autocorrelation function within subspaces followed an inverse exponential, leading to a *Toeplitz* structure in the subspace correlation matrix. The maximal correlation in a subspace varied from 0 to 0.65, and was different for each subspace in the same experimental run. In each of the ten runs for each of the six correlation cases considered, a new, unique square mixing matrix A (`V = C`) with condition number `cond(A) = 3` was randomly generated. Likewise, each run was initialized with a different random row-orthogonal `W_0` (60 runs total). Finally, we compared the performance of IVA-GL (i.e., IVA-G with Newton step (see `icatb_iva_second_order.m` in `other_methods` folder) as an initializer for IVA-L  with relative gradient (see `icatb_iva_laplace.m` in `other_methods` folder)) and MISA in each setting.

`Optimization parameters`: No data reduction was employed. For MISA, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of `W` were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 10000, the maximum number of
iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for `L-BFGS-B` was set to 1. Also, `TolX = 10^-30` and `TolFun = 10^-4`. For IVA-G, we utilized all the defaults in [GIFT] and `opt approach` was set to `newton`. For IVA-L, we utilized all the defaults in [GIFT], `whiten` was set to `false`, `maxIter` was set to 1024, and `InitialW` was set the unmixing matrix estimated with IVA-G.

### `IVA2` (`V < N` case):
Here we assess the effect of additive sensor noise and condition number in a larger IVA problem with rectangular mixing matrix `A` and an abundant number of observations `N`. The experiment was setup with `M = 16` datasets, each containing `C_m = 75` sources for a total of `K = C_m = 75` sixteen-dimensional subspaces, and `N = 66000` observations sampled
independently from a multivariate Laplace distribution. Similarly to experiment `IVA1`, an inverse exponential autocorrelation function with maximal correlation varying from 0 to 0.5 was selected. In each of the ten instances of this experiment, a new, unique (`V x C`) *rectangular* mixing matrix `A` (`V = 250`) was randomly generated. For each instance, ten runs were performed, each with a different random row-orthogonal `W_0` initialization (thus, 100 runs per experimental condition). The
experiments were broken into two parts, a) and b), exactly as described in experiment `ICA1`. We also evaluate the performance of data-reduction by group PCA (gPCA) or PRE, followed by either IVA-L or MISA.

`Optimization parameters`: For the data reduction, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of `W` were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 50000, the maximum number of iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for L-BFGS-B was set to 5. `TolX` and `TolFun` were set as described in Section II-B of the supplemental material for the acompanying publication, with `b = 80` for part a), and b = 4 for part b). For MISA, we utilized the same parameters, except: `TolFun = 10^-4`, `TolX = 10^-9`, the initial barrier parameter was set to 1, and the number of past gradients utilized for L-BFGS-B was set to 10. For IVA-L, we utilized all the defaults in [GIFT] and `whiten` was set to `false`, `alpha0` was set to 1, `terminationCriterion` was set to `ChangeInW`, and `maxIter` was set to `min(10000, 4*MISAiter)`, where
`MISAiter` is the number of iterations until convergence for MISA on the same problem, from the same starting point.

### `ISA1_2` (`V < N`, `V = C` case):
In these experiments, we want to assess the performance in independent subspace analysis (ISA) problems when no data reduction is required (i.e., `V = C`), noise is absent, and the number of observations `N` is abundant. The `ISA1` experiment was setup with `M = 1` dataset containing `C = C_1 = 28` sources organized into `K = 7` `dk`-dimensional subspaces (see `dk` values below), with `N = 32968` observations sampled independently from a multivariate Laplace distribution. Correlation was absent within all subspaces (i.e., only higher-order statistics (HOS)
dependence was present). In each of the ten runs for each of the two cases below, a new, unique square mixing matrix `A` (`V = C`) with condition number `cond(A) = 3` was randomly generated. Likewise, each run was initialized with a different random row-orthogonal `W_0` (20 runs total). The experiments were broken into two cases:

- Case 1: subspace sizes `dk = k`, `k = [1, 2, 3, 4, 5, 6, 7]`.
- Case 2: subspace sizes `dk = 4` for all subspaces.

The `ISA2` experiment was identical to experiment `ISA1`, except with a within-subspace correlation structure like the one described in experiment `IVA1`, with maximal correlation varying from 0.2 to 0.75. We also evaluate the performance of [JBD-SOS], [isa_est] (in Case 2 only) and MISA-GP in each setting.

`Optimization parameters`: No data reduction was employed. For MISA-GP, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of W were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 10000, the maximum number of
iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for L-BFGS-B was set to 1. Also, TolX = 10^-30 and TolFun = 10^-4. For [JBD-SOS], we provided 500 data covariance matrices as input (computed after PRE data reduction), `threshold` was set to 10^-9, and `max sweep` was set to 512.

[JBD-SOS]: https://www.irit.fr/~Dana.Lahat/jbd.zip
[isa_est]: http://ai.stanford.edu/~quocle/video_release.tar.gz

### `ISA3` (`V > N` case):
Here we assess the effect of additive sensor noise and condition number in a mildly large ISA problem with rectangular mixing matrix `A` and a fairly low number of observations `N`. The experiment was setup with `M = 1` dataset containing `C = C_1 = 51` sources organized into `K = 18` `dk`-dimensional subspaces, with `dk = [1:5, 5:1, 1:5, 2, 2, 2]` (where
`1:5` means “one through 5”) and `N = 5250` observations sampled independently from a multivariate Laplace distribution. Similarly to experiment `IVA1`, an inverse exponential autocorrelation function with maximal correlation varying from 0 to 0:5 was selected. In each of the ten instances of this experiment, a new, unique (`V x C`) *rectangular* mixing matrix `A` (`V = 8000`) was randomly generated. For each instance, ten runs were performed, each with a different random row-orthogonal `W_0`
initialization (thus, 100 runs per experimental condition). The experiments were broken into two parts, a) and b), exactly as described in experiment `ICA1`. We also evaluate the performance of [JBD-SOS], MISA, and MISA-GP in each setting.

`Optimization parameters`: For the data reduction, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of `W` were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 50000, the maximum number of iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for L-BFGS-B was set to 5. `TolX` and `TolFun` were set as described in Section II-B of the supplemental material for the acompanying publication, with `b = 80` for part a), and `b = 4` for part b). For MISA-GP, we utilized the same parameters, except: `TolFun = 10^-4`, `TolX = 10^-9`, the initial barrier parameter was set to 1, and the number of past gradients utilized for L-BFGS-B was set to 15. For [JBD-SOS], we provided 500 data covariance matrices as input (computed after PRE data reduction), `threshold` was set to 10^-6, and `max sweep` was set to 512.

## Hybrid data examples
### `simhybridICAfMRI.m`

This example considers *single-subject temporal ICA of fast acquisition fMRI*. Details are outlined in the accompanying publication.

`Optimization parameters`: For the data reduction with PRE, estimating `A` from `pinv(W)`, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of `W` were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 50000, the maximum number of iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for L-BFGS-B was set to 5. `TolX` and `TolFun` were set as described in Section II-B of the supplemental material for the acompanying publication, with `b = 4`, or as 10^-8 (whichever was smallest). This optimization was performed on each dataset separately. For MISA, we utilized the same parameters, except: `TolFun = TolX = 10^-9`, the initial barrier parameter was set to 1, and the number of past gradients utilized for L-BFGS-B was set to 10.

### `simhybridMMIVA.m`

This example considers *multimodal IVA of structural MRI, functional MRI, and Fractional Anisotropy (FA) difusion MRI*. Details are outlined in the accompanying publication.

`Optimization parameters`: For the data reduction with PRE, estimating `A` from `pinv(W)`, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of `W` were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 50000, the maximum number of iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for L-BFGS-B was set to 5. `TolX` and `TolFun` were set as described in Section II-B of the supplemental material for the acompanying publication, with `b = 4`, or as 10^-8 (whichever was smallest). This optimization was performed on each dataset separately. For MISA, we utilized the same parameters, except: `TolFun = TolX = 10^-9`, the initial barrier parameter was set to 1, and the number of past gradients utilized for L-BFGS-B was set to 10.


### `simhybridMISA.m` 
This example illustrates the value of *MISA without data reduction* for *multimodal fusion of electroencephalography (EEG) event-related potentials (ERP) and fMRI* datasets. Details are outlined in the accompanying publication.

`Optimization parameters`: For the data reduction with simple RE, estimating `A` from `W^T`, we utilized MATLAB’s Optimization Toolbox function `fmincon.m`. The bounds on the values taken by each element of `W` were set to -100 and 100, their typical value was set to `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 50000, the maximum number of iterations was set to 10000, the initial barrier parameter was set to 0.1, and the number of past gradients utilized for L-BFGS-B was set to 5. `TolX` and `TolFun` were set as described in Section II-B of the supplemental material for the acompanying publication, with `b = 4`. This optimization was performed on each dataset separately.
For MISA-GP, we utilized the same parameters, except: `TypicalX = 0.1`, the maximal number of objective function evaluations was set to 50000, the maximum number of iterations was set to 150, the initial barrier parameter was set to 1, the number of past gradients utilized for L-BFGS-B was set to 1 (which is equivalent to conjugate gradient), and `TolFun = TolX = 0.5*N*10^-9`.
For the following corrective ICAs, we also used MISA with the latter set of parameters.

# Supporting material

## Folders
- `toy_example`: data and scripts for the toy example in the main page.
- `MCIv4`: a suite of functions for nice display of spatial maps in MATLAB. Written by [Elena Allen].

[Elena Allen]: https://scholar.google.com/citations?user=vrXdBU4AAAAJ&hl=en

## Files
- `setup_*_MISAKRE.m`: simplified interface for instantiation of MISAKRE objects using simulated or hybrid data.
- `sim_basic_*.m`: simplified interface for generating simulated data.
- `sim_hybrid_*.m`: simplified interface for generating hybrid data.
- `doSecondPCAStep.m`: supporting function for group PCA.
- `my_3views.m`: supporting function for display.
- `to_vol.m`: supporting function for display.
