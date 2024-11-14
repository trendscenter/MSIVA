import os
import mat73
import numpy as np
import scipy.io as sio
from scipy import stats
import statsmodels
import statsmodels.api as sm

def mask(arr):
    ind = ~np.isnan(arr)
    arr_ = arr[ind]
    return arr_, ind

dirpath = "/data/users4/xli/MSIVA/mat"
ver = 4
n_var = 25

phenotype_path = os.path.join(dirpath, f"ukb_2907sub_{n_var}var.mat") 
phenotype = sio.loadmat(phenotype_path)['phenotype_value']

delta2p = mat73.loadmat( os.path.join(dirpath, f"v{ver}", "delta2p_raw", "delta2p_sMRI.mat") )["delta2p"]
n_subject = delta2p.shape[0]
n_predictor = delta2p.shape[1]
n_voxel = delta2p.shape[2]
n_phenotype = phenotype.shape[1]

rall_list, pall_list = [], []
for i in range(n_phenotype):
    r = np.zeros((n_voxel, n_predictor))
    p = np.zeros((n_voxel, n_predictor))
    pheno, ind = mask(phenotype[:,i])
    for v in range(n_voxel):
        for scv in range(n_predictor):
            delta = delta2p[:,scv,v][ind]
            r[v,scv], p[v,scv] = stats.pearsonr(delta, pheno)
    sio.savemat(os.path.join(dirpath, f"v{ver}", f"phenotype_map_{n_var}var", f"phenotype_map{i+1}.mat"), {"r":r, "p":p})
    rall_list.append(r)
    pall_list.append(p)

# FDR correction
pall = np.dstack(pall_list)
num_vox = pall.shape[0]
num_scv = pall.shape[1]
num_var = pall.shape[2]
print(num_vox, num_scv, num_var)

rejected_1d, p_corrected_1d = statsmodels.stats.multitest.fdrcorrection(pall.flatten(), alpha=0.05, method='poscorr', is_sorted=False)
p_corrected = p_corrected_1d.reshape( (num_vox, num_scv, num_var) )
rejected = rejected_1d.reshape( (num_vox, num_scv, num_var) )

rall = np.dstack(rall_list)
sio.savemat(os.path.join(dirpath, f"v{ver}", f"phenotype_map_{n_var}var", "pearsonr.mat"), {"r": rall})
sio.savemat(os.path.join(dirpath, f"v{ver}", f"phenotype_map_{n_var}var", "p_corrected.mat"), {"p_corrected": p_corrected})
sio.savemat(os.path.join(dirpath, f"v{ver}", f"phenotype_map_{n_var}var", "rejected.mat"), {"rejected": rejected})