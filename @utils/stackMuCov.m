function w = stackMuCov(mu_, cov_)

w = mu_;
for ii = 1:size(cov_,1)
    w = [w; cov_(ii,ii:size(cov_,1))'];
end

end