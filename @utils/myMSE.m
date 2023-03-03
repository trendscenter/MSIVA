function MSE = myMSE(O,Rr)

d = min(size(Rr));
% MSE = 1 - (1/d)*trace(Rr*Rr') + O.myMD(Rr)^2;
Rr = abs(Rr);
D = zeros(d,d);
for ii = 1:d
    for jj = 1:d
        tmp = - Rr(ii,:);
        tmp(jj) = tmp(jj) + 1;
        D(ii,jj) = sum(tmp.^2);
    end
end
p = O.munkres(D);
Rr = Rr(:,p);
% figure, imagesc(Rr,[-1 1])

MSE = 2 - (2/d)*trace(Rr);

end

% function MD = myMD(WAr)
% 
% d = min(size(WAr));
% G = WAr.^2;
% p = munkres(bsxfun(@times,G,1./sum(G,2)));
% WAr = WAr(:,p);
% 
% MD = sqrt(d - trace(WAr)/(d-1));
% 
% end