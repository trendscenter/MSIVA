function [out, WAr] = MMD(O,W,A,S)

WA = cellfun(@(w,a) w*a, W, A, 'Un', 0);
WA = blkdiag(WA{:});

K = size(S{1},1);
Stmp = [S{:}];
idx = [];
cps = [];
for kk = 1:K
    tidx = find(Stmp(kk,:));
    idx = [idx tidx];
    cps(kk) = length(tidx);
end

WA = WA(idx,idx);

% cps = sum(Stmp,2); % This is bugged... crashes the server for large number of datasets...

% figure, imagesc(WA)
WAr = mat2cell(WA,cps,cps);
%WAr = cell2mat(cellfun(@(wa) mean(abs(wa(:))), WAr, 'Un', 0));
WAr = cell2mat(cellfun(@(wa) sum(abs(wa(:))), WAr, 'Un', 0));

out = O.myMD(WAr);

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