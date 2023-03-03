function [mu_, cov_] = unstackMuCov(w,M,V)

% Unstack wout: *wout is a col. vector. *Number of elements in Wf = Cf * Nf, different for each f
numel_Wf = [V(M), V(M)*(V(M)+1)/2];
% To partition wout into M Wf's, need start/end indexes corresponding to each f.
ixe = cumsum(numel_Wf);       % End index
ixb = 1 + [0 ixe(1:(end-1))]; % Begin index
% Parse wout:  Elements of wout are rows of Wf (M=1 first).
mu_ = w(ixb(1):ixe(1));
cov_long = w(ixb(2):ixe(2));
cov_ = zeros(V(M));
ll = 0;
for ii = 1:V(M)
    for jj = ii:V(M)
        ll = ll + 1;
        cov_(ii,jj) = cov_long(ll);
        if ii ~= jj
            cov_(jj,ii) = cov_long(ll);
        end
    end
end