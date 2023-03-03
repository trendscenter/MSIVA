function [out, Rr] = MMSE(O,Y,Yh,S)

Ryyh = cellfun(@(yh,y) corr(y',yh'), Y, Yh, 'Un', 0);
Ryyh = blkdiag(Ryyh{:});

K = size(S{1},1);
Stmp = [S{:}];
idx = [];
cps = [];
for kk = 1:K
    tidx = find(Stmp(kk,:));
    idx = [idx tidx];
    cps(kk) = length(tidx);
end

Ryyh = Ryyh(idx,idx);

% cps = sum(Stmp,2); % This is bugged... crashes the server for large number of datasets...

% figure, imagesc(WA)
Rr = mat2cell(Ryyh,cps,cps);
%WAr = cell2mat(cellfun(@(wa) mean(abs(wa(:))), WAr, 'Un', 0));
Rr = cell2mat(cellfun(@(r) max(abs(r(:))), Rr, 'Un', 0));

% out = O.myMSE(bsxfun(@rdivide,Rr,max(Rr,[],2)));
out = O.myMSE(Rr);

end

% function MSE = myMSE(Rr)
% 
% d = min(size(Rr));
% MSE = 1 - (1/d)*trace(R*R') + myMD(Rr)^2;
% 
% end
% 
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