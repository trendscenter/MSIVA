function [Y1wht,whtM,dewhtM,Y1br,brM,debrM] = myPCA(Y1)

V = size(Y1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA of Y1, with ordering by descending variance:
% Covariance matrix of Y1
C1 = (Y1*Y1')./(V-1);
% SVD Decomposition of C1
%[U,s,V] = svd(C1);
% Eigval Decomp of C1
[V1,E1] = eig(C1);
% Reorder cols of E1 and V1 (largest to smallest)
[temp, seq] = sort(diag(E1),'descend');
E1 = E1(seq,seq);
V1 = V1(:,seq);
% Rotate: project Y1 onto eig-vector space V1
Y1rot = V1'*Y1;
% Whitening: rescale Y1rot for unit variances (in addition to rotation)
Y1wht = sqrt(E1)\Y1rot; % Same as: Ywht = inv(sqrt(E))*Yrot;
whtM = sqrt(E1)\(V1');
dewhtM = V1*sqrt(E1);
% Back-recontruction: rotate back to data space
Y1br = V1*Y1wht;
brM = V1;
debrM = V1';
