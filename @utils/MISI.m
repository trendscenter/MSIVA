function [out, WAr] = MISI(O,W,A,S)

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

out = O.myISI(WAr);

end

% function ISI = myISI(WA)
% 
% N = length(WA);
% 
% WA = abs(WA);
% 
% ISI = 0;
% ISI = ISI + sum(sum(WA,2)./max(WA,[],2) - 1);
% % max(WA,[],2)
% ISI = ISI + sum(sum(WA,1)./max(WA,[],1) - 1);
% % max(WA,[],1)
% 
% ISI = ISI/(2*N*(N-1));
% 
% end