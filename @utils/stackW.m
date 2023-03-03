function w = stackW(W)

w = cellfun(@(x) reshape(x',1,[]),W,'Un',0);
w = [w{:}]';

end