function MD = myMD(O,WAr)

d = min(size(WAr));
G = WAr.^2;
G = bsxfun(@times,G,1./sum(G,2));
D = zeros(d,d);
for ii = 1:d
    for jj = 1:d
        tmp = - G(ii,:);
        tmp(jj) = tmp(jj) + 1;
        D(ii,jj) = sum(tmp.^2);
    end
end
p = O.munkres(D);
G = G(:,p);

MD = sqrt((d - trace(G))/(d-1));

end