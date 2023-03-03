function update(O,S,M,b,l,e)

O.S = S;
O.M = M;

O.d = full(sum([O.S{O.M}],2));  % Dimensionality of each subspace
O.nes = O.d~=0;             % Non-empty subspace indexes
O.d_k = cellfun(@(s) full(sum(s,2)), O.S,'Un',0);

O.K = sum(O.nes);

rs = ones(size(O.S{O.M(1)},1), 1); % rows of subspace

if length(b) == 1
    O.beta = b*rs;
else
    O.beta = b;
end

if length(e) == 1
    O.eta = e*rs;
else
    O.eta = e;
end

O.nu = (2*O.eta + O.d - 2)./(2*O.beta);

O.auto_tune('lambda', l);
% O.lambda = l;
O.a = (O.lambda.^(-1./(O.beta)) .* gamma(O.nu + 1./O.beta)) ./ ...
    (O.d .* gamma(O.nu));

end