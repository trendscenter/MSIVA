function sim1 = sim_basic_IVA(seed,K,V,M_Tot,N,Acond,SNR)
% MISA
% seed    = [];
%N       = 200000;
% % S_ = {};
% for kk = 1:length(C)
%     S_{kk} = sum(C(1:(kk-1)))+(1:C(kk));
% end
% S_ = S_';
% [K, M_Tot] = size(S_);
S_ = num2cell(repmat((1:K)', 1, M_Tot));
C = [];
for mm = 1:M_Tot
    C = [C sum([S_{:,mm}] ~= 0)];
end
cr = linspace(.65,.85,K);
for kk = 1:K
    dist_params(kk).name = 'mvn';
    dist_params(kk).mu   = zeros(sum([S_{kk,:}] ~= 0),1);
    dist_params(kk).CORR = cr(kk); % Ignored for ICA
end
%V       = C;
Atype   = 'Generated';
%Acond   = 3;
A       = {};
%SNR     = (1+99)/99; % (data power + noise power) / noise power 
predictors = cell(1,M_Tot);
regtype = 'LS';
sim1    = gsd(seed, N, M_Tot, K, C, S_, dist_params, ...
    V, Atype, Acond, A, SNR, predictors, regtype);
