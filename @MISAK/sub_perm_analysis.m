function [w0, shuff] = sub_perm_analysis(O, w0)

M = O.M; % Save current M
S = O.S; % Save current network structure S
% % Save parameters
% old_beta = O.beta;
% old_lambda = O.lambda;
% old_eta = O.eta;
% old_d_k = O.d_k;
O.updatesc(false);

% Subspace permutation analysis:
% List all possible subspace permutations, based on original (user supplied) S
subspace_perms = gen_list_Gperms(S,M);
spval = [];                                       % Container for cost value
for kk = 1:size(subspace_perms{M(1)},1)           % For each subspace perm
    % Apply perm to each dataset
    S_ = cellfun(@(s,p) s(p(kk,:),:), S(M), subspace_perms(M), 'Un', 0);
    O.update(S_,M,O.beta,O.lambda,O.eta);
    spval(kk) = O.objective(w0);                  % Compute cost
end
[~, the_min] = min(spval);                        % best subspace permutation
S_ = cellfun(@(s,p) s(p(the_min,:),:), S(M), subspace_perms(M), 'Un', 0);
% figure, imagesc(full(S_{1}),[0 1]), axis equal tight; title('AFTER subspace match: ESTIMATED')
% After subspace permutation, components are no longer arranged in
% the same order as subspaces. This will reshuffle the components.
shuff = cell(1,max(M));
for mm = M
    shuff{mm} = [];
    for ss = 1:size(S_{M(1)},1)                   % Loop through suspaces
        shuff{mm} = [shuff{mm} find(S_{mm}(ss,:))];
    end
end

% Apply shuffle to original structure
O.updatesc(true);
O.update(S,M,O.beta,O.lambda,O.eta);
O.objective(O.ut.stackW(cellfun(@(w,s) w(s,:), O.W(M), shuff(M), 'Un', 0)));
w0 = O.ut.stackW(O.W);

    function out = gen_list_Gperms(Sdes,M_)
        
        mygperms = cell(1,max(M_));
        for ff = M_
            Sd = full(sum(Sdes{ff},2));             % Number of components in each subspace
            
            [t1, t2, t3] = unique(Sd, 'stable');    % Unique subspace sizes
            stot = zeros(size(t1));
            for tt = 1:length(t1)
                % Total number of components in subspaces of size t1
                stot(tt) = sum(Sd == t1(tt));
            end
            
            ix = cell(length(t1),1);
            for xx = 1:length(ix)
                % All possible permutations of subspaces with same size
                ix{xx} = perms(find(Sd == t1(xx))');
            end
            % Compile together all subspace permutations from all subspaces
            seqs = myGperms(ix);
            
            %
            tokeep = 1:length(Sd);                  %
            shuf = tokeep;                          % Current subspace shuffle (no shuffle)
            se = 0;                                 % end index
            cumstot = cumsum(stot);                 % Total number fo components
            for kk = 1:length(t1)                   % For each subspace size
                sb = se + 1;                        % beginning index
                se = cumstot(kk);                   % end index: select components from sb to se
                shuf(t3 == kk) = tokeep(sb:se);     % subspasce indices in order of size t1
            end
            % Shuffle back into natural occurrence order
            mygperms{ff} = seqs(:,shuf);
            
        end
        
        % Compile together all subspace permutations from all datasets
        mygperms = myGperms(mygperms(M_));          % datasets stacked side by side, rows are perms
        
        K = size(Sdes{M_(1)},1);
        
        ixb = 1:K:(K*length(M_));                   % Begin indexes
        ixe = K:K:(K*length(M_));                   % End indexes
        
        out = cell(1,max(M_));
        for kk = 1:length(M_)
            % split mygperms into M_ blocks, 1 per dataset
            % out contains lists of permutations for each dataset
            out{M_(kk)} = mygperms(:,ixb(kk):ixe(kk));
        end
        
    end

    function out_ = myGperms(list)
        % Compile together all subspace permutations from all subspaces
        % list contains permutations for each subspace
        out_ = [];
        if length(list) > 1
            for kk = 1:size(list{1},1)
                % Combine items of lists
                out_ = [out_; permme(list{1}(kk,:),list(2:end))];
            end
        else
            out_ = list{1};
        end
    end

    function out__ = permme(L, list)
        % Combine items of lists
        out__ = [];
        if isempty(list)
            out__ = L;
        else
            for kk_ = 1:size(list{1},1)
                % Combine items of lists
                out__ = [out__; permme([L list{1}(kk_,:)],list(2:end))];
            end
        end
        
    end

end