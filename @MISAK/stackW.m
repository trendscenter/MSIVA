function w = stackW(O,W)

% w = cellfun(@(x) reshape(x',1,[]),W,'Un',0);
% w = [w{:}]';
w = O.ut.stackW(W);

end