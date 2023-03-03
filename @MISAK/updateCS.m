function updateCS(O,S)

for mm = O.M
    new_cols = size(S{mm},2);
    old_cols = size(O.S{mm},2);
    if new_cols > old_cols
        
        A = O.W{mm}';
        [Q,~] = qr(A,0);
        B = [Q randn(size(A,1), new_cols - old_cols)];
        [Q_,~] = qr(B,0);
        Q_ = Q_(:,old_cols+1:end);
        
        tmp = zeros(size(B'));
        for kk = find(O.nes)'
            
            old_cc = find(O.S{mm}(kk,:) == 1);
            new_cc = find(S{mm}(kk,:) == 1);
            df = length(new_cc) - length(old_cc);
            
            if df <= 0
                tmp(new_cc,:) = A(:,old_cc(1:length(new_cc)))';
            else
                aux = B(:,old_cc)';
                mx = mean(aux);
                new_b = Q_(:,1:df)';
                Q_ = Q_(:,(df+1):end);
                aux = [aux; bsxfun(@plus,(.5*max(abs(mx))/max(abs(new_b(:)))) * new_b, mx)];
                tmp(new_cc,:) = aux;
            end
            
        end
        O.W{mm} = tmp;
        
    else
        
        ixs = [];
        for kk = find(O.nes)'
            ixs = [ixs find(O.S{mm}(kk,:) == 1, 1)];
        end
        O.W{mm} = O.W{mm}(ixs,:);
        
    end
    
end
O.S = S;
O.C = cellfun(@(s) size(s,2), O.S); % Number of components per dataset

% Compute source estimates (network output):
O.Y = cell(1,max(O.M));
O.Y(O.M) = cellfun(@mtimes, O.W(O.M), O.X(O.M), 'Un', 0);

O.d = full(sum([O.S{O.M}],2));  % Dimensionality of each subspace
O.nes = O.d~=0;             % Non-empty subspace indexes
O.d_k = cellfun(@(s) full(sum(s,2)), O.S,'Un',0);

O.K = sum(O.nes);

end