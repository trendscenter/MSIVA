import mat73
import numpy as np
import scipy.io as sio
import hdmedians as hd

ver = 4
delta2p = mat73.loadmat(f"/data/users4/xli/MSIVA/mat/v{ver}/delta2/delta2p_sMRI.mat")["delta2p"] # subject x predictor x voxel

num_predictor = delta2p.shape[1]
num_voxel = delta2p.shape[2]

geomedian = np.zeros((num_predictor, num_voxel))
for p in range(num_predictor):
    geomedian[p,:] = hd.geomedian(delta2p[:,p,:], axis=0)

sio.savemat(f"/data/users4/xli/MSIVA/mat/v{ver}/delta2p/delta2p_geomedian_sMRI.mat", {"geomedian": geomedian})