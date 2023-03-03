function sim1 = sim_hybrid_IVA(seed,K,N,A,SNR,distname)

M_Tot = length(A);
V = cellfun(@(x) size(x,1), A);
S_ = num2cell(repmat((1:K)', 1, M_Tot));
C = [];
for mm = 1:M_Tot
    C = [C sum([S_{:,mm}] ~= 0)];
end
cr = linspace(.65,.85,K);
for kk = 1:K
    dist_params(kk).name = distname;
    dist_params(kk).mu   = zeros(sum([S_{kk,:}] ~= 0),1);
    dist_params(kk).CORR = cr(kk); % Ignored for ICA
end
Atype   = 'Provided';
Acond   = [];
predictors = cell(1,M_Tot);
regtype = 'LS';
sim1    = gsd(seed, N, M_Tot, K, C, S_, dist_params, ...
    V, Atype, Acond, A, SNR, predictors, regtype);
