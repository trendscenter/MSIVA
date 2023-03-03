function [wM, groupPCASig] = doSecondPCAStep(data,numComp)
%% Do second step pca for input to IVA

% numComp = size(data, 1);
numVoxels = size(data, 2);
numSubjects = size(data, 3);

data = permute(data, [2, 1, 3]);
data = reshape(data, size(data, 1), size(data, 2)*size(data, 3));
[groupPCASig, dewhiteM] = icatb_calculate_pca(data, numComp, 'type', 'mpowit', 'whiten', true);
wM = pinv(dewhiteM);