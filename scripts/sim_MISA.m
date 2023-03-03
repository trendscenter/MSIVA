function sim1 = sim_MISA(seed,S_,V,N,Acond,SNR)
% sim_MISA(seed,K,V,M_Tot,N,Acond,SNR)
% K: number of subspaces -> a vector of number of sources per subspace
% V: data dimension
% M_Tot: total number of modalities
% N: number of samples
% Acond: condition number of mixing matrix (ratio of largest to smallest singular value in the matrix), 
% if it's too large, can't invert the matrix
% SNR: signal to noise ratio

% MISA
% seed    = [];
%N       = 200000;
% % S_ = {};
% for kk = 1:length(C)
%     S_{kk} = sum(C(1:(kk-1)))+(1:C(kk));
% end
% S_ = S_';
% [K, M_Tot] = size(S_);

% S_ = mat2cell(repmat((1:sum(K)), M_Tot, 1), ones(1,M_Tot), K)';

[lenK, M_Tot] = size(S_);

% C: vector of total number of sources per modality, infer from S_
C = [];
for mm = 1:M_Tot
    C = [C sum([S_{:,mm}] ~= 0)];
end

% cr: maximum value of the correlation for each subspace
% MVK: a mixture of Gaussian where weights follow a gamma distribution
cr = linspace(.85,.65,lenK);
for kk = 1:lenK
    dist_params(kk).name = 'mvl';
    dist_params(kk).mu   = zeros(sum([S_{kk,:}] ~= 0),1);
    dist_params(kk).CORR = cr(kk); % Ignored for ICA
end

% V = C; 
Atype = 'Generated';
% Acond = 3;
A = {};
% SNR = (1+99)/99; % (data power + noise power) / noise power 
predictors = cell(1,M_Tot);
regtype = 'LS'; % not used
sim1 = gsd(seed, N, M_Tot, lenK, C, S_, dist_params, ...
    V, Atype, Acond, A, SNR, predictors, regtype);

