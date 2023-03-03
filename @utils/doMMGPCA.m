function [whtM, H] = doMMGPCA(X, comps, rec_type)

M = 1:length(X);
N = size(X{1},2);

V = cellfun(@(x) size(x,1), X, 'Un', 0);

if N <= min([V{:}])
    % Xm.T * Xm, weighted avg per modality --> cov(X)
    % Equivalent to X_cat.T*X_cat, with some scaling --> cov(X_cat)
    % U,S,H = SVD(X_cat), H,S2 = EIG(X_cat.T*X_cat)
    % Then S2 = S.^2
    
    cvx = zeros(N);
    
    for mm = M
        cvx_ = cov(X{mm});
        cvx = cvx + cvx_./(length(M)*trace(cvx_)/N);
    end
    
    % Subject-level PCA reduction...
    [H,lambda] = eigs(cvx,comps);

else
    % U,S,H = SVD(X_cat), U,S2 = EIG(X_cat*X_cat.T)
    % scale data and concatenate
    X_cat = cell2mat(cellfun(@(x) sqrt(N/(length(M)*sum(x(:).^2))) * x, X', 'Un', 0));
    cvx = X_cat*X_cat';

    [U,lambda] = eigs(cvx,comps);
    H = ((diag(1./sqrt(diag(lambda))) * U') * X_cat)';
end


A = cellfun(@(x) sqrt(N./(length(M)*sum(x(:).^2)))*(x*H), X, 'Un', 0);
norm_A = cellfun(@(a) sum(a.^2), A, 'Un', 0)';
if M == 1
    norm_A = sqrt(cell2mat(norm_A));
else
    norm_A = sqrt(sum(cell2mat(norm_A)));
end
A = cellfun(@(a) a./repmat(norm_A,size(a,1),1),A,'Un',0);

if strcmpi(rec_type, 'WT')
    whtM = cellfun(@(a) a', A, 'Un', 0);
elseif strcmpi(rec_type, 'PINV')
    whtM = cellfun(@(a) pinv(a), A, 'Un', 0);
elseif strcmpi(rec_type, 'REG')
    % To Do...
end

whtM = cellfun(@(x,w) sqrt(N-1) * sqrt(N./(length(M)*sum(x(:).^2))) * repmat(1./(norm_A'),1,size(w,2)) .* w, X, whtM, 'Un', 0);
% whtM = cellfun(@(x,w) sqrt(N-1) * sqrt(N./(length(M)*sum(x(:).^2))) * diag(1./(norm_A')) * w, X, whtM, 'Un', 0);

H = sqrt(N-1) * H';