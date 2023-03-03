function sim1 = sim_basic_ICA(seed,C,V,N,Acond,SNR)
% ICA
% seed    = [];
%N       = 200000;
S_      = mat2cell((1:C)', ones(C,1));
[K, M_Tot] = size(S_);
C = [];
for mm = 1:M_Tot
    C = [C sum([S_{:,mm}] ~= 0)];
end
for kk = 1:K
    dist_params(kk).name = 'mvl';
    dist_params(kk).mu   = zeros(sum([S_{kk,:}] ~= 0),1);
    dist_params(kk).CORR = 0.5; % Ignored for ICA
end
%V       = C;
Atype   = 'Generated';
%Acond   = 3;
A       = {};
%SNR     = (1+99)/99; % (data power + noise power) / noise power 
predictors = 3:4;%{1:2, [], 1:2};
regtype = 'LS';
sim1    = gsd(seed, N, M_Tot, K, C, S_, dist_params, ...
    V, Atype, Acond, A, SNR, predictors, regtype);
