function w0 = greedysearch(O,varargin)

if ~isempty(varargin)
    sim1 = varargin{1};
end

if length(O.M) ~= 1 || ... % More than 1 dataset
        O.K ~= O.C(O.M(1)) || ... % If 1 dataset, subspace size K ~= num comp C
        prod(diag((O.S{O.M(1)}(O.nes,1:O.C(O.M(1))))) == ones(O.C(O.M(1)),1)) ~= 1
    % If 1 dataset, and K == C, S is not identity mtx
    
    % Run ICA on each dataset separately (assuming MISA already ran
    % once):
    M = O.M; % Save current M
    S = O.S; % Save current network structure S
    S_ = S;
%     REl = O.RElambda;
    % Save parameters
    old_beta = O.beta;
    old_lambda = O.lambda;
    old_eta = O.eta;
    old_d_k = O.d_k;
    sc = O.sc; % Save scale-control parameter
    oldW = O.W; % Save surrent W
    for mm = M
        % Set S to sparse identity matrix
        S_{mm} = sparse(1:O.C(mm),1:O.C(mm),ones(O.C(mm),1),O.C(mm),O.C(mm),O.C(mm)); %eye(O.C(mm));
        % redefine parameters: accounts for changing subspace structure/size
        ixe = cumsum(cell2mat(old_d_k(mm)));
        ixb = [1; (ixe(1:end-1)+1)];
        b = []; l = []; e = [];
        for kk = 1:length(ixe)
            b(ixb(kk):ixe(kk)) = old_beta(1);
            l(ixb(kk):ixe(kk)) = old_lambda(1);
            e(ixb(kk):ixe(kk)) = old_eta(1);
        end
        O.update(S_,mm,b',l',e'); % Update S and set M = mm
        % Initial w for ICA is the one resulting from the previous MISA
        % run:
        w0 = O.ut.stackW(O.W(mm));
        
        % Setup and run ICA on dataset mm:
%         O.updateRElambda(0);
%         O.updateRElambda(O.opt_RE(w0));
        optprob = O.ut.getop(w0, @(x) O.objective(x), [], 1, {'lbfgs' 10}, 1e-9); % For sim results, remove barrier, lbfgs and tol params
        optprob.options.MaxIter = 1000; % For sim results, comment this
%         optprob.options.Display = 'iter-detailed';
%         optprob.options.PlotFcns = {@optimplotconstrviolation @optimplotstepsize @optimplotfirstorderopt @optimplotfval};
        [woutT,fvalT,exitflagT,outputT] = fmincon(optprob);
        %W = O.ut.unstackW(woutT,O.M,O.C,O.V);
        %WoutT{mm} = W{mm};
    end
    
    % Loop through all components in all datasets:
    O.updatesc(false);                                % Turn off scale-control
    for mm = M
        w0 = O.ut.stackW(O.W(mm));
        ixe = cumsum(cell2mat(old_d_k(mm)));
        ixb = [1; (ixe(1:end-1)+1)];
        b = []; l = []; e = [];
        for kk = 1:length(ixe)
            b(ixb(kk):ixe(kk)) = old_beta(1);
            l(ixb(kk):ixe(kk)) = old_lambda(1);
            e(ixb(kk):ixe(kk)) = old_eta(1);
        end
        O.update(S_,mm,b',l',e');
        
        % Component permutation analysis (per dataset):
        % The following is a GREEDY SEARCH on the space component
        % assignements to subspaces
%         figure
        for cc = randperm(size(S{mm}, 2))%1:size(S_{mm}, 2)%
            current = find(S_{mm}(:,cc));             % current subspace for component cc
            kk = find(S_{mm}(current,:));
            S_{mm}(:,kk) = zeros(size(S_{mm},1),length(kk));   % empty column cc
            misa_values = [];                         % container for cost values
            for ss = 1:size(S_{mm},1)
                if ss ~= 1
                    S_{mm}(ss-1,kk) = 0;              % Remove previous assignment
                end
                S_{mm}(ss,kk) = 1;                    % Assign component cc to subspace ss
                O.update(S_,mm,b',l',e');
                misa_values(ss) = O.objective(w0);    % Compute cost
            end
            S_{mm}(ss,kk) = 0;                        % Remove previous assignment
            S_{mm} = [S_{mm}; zeros(1,size(S_{mm},2))]; % Try adding 1 more subspace
            S_{mm}(end,kk) = 1;                       % Assign to new subspace
            O.update(S_,mm,[b'; b(end)],[l'; l(end)],[e'; e(end)]); % Update S and set M = mm
            misa_values(ss+1) = O.objective(w0);
            
%             global sim_misa;
%             global W;
%             figure,imagesc(O.W{mm}*W{mm}*sim_misa.A{mm},max(max(abs(O.W{mm}*W{mm}*sim_misa.A{mm}))).*[-1 1]);colorbar();
%             figure,imagesc(O.W{end}*W{end}*sim_misa.A{end},max(max(abs(O.W{end}*W{end}*sim_misa.A{end}))).*[-1 1]);colorbar();

            [~,ix] = min(misa_values);                % best subspace for component cc
%             imagesc(full(S_{mm}))
            drawnow()
            ix = find(misa_values == misa_values(ix));
            if sum(ix == (ss+1)) > 0
                ix = ss+1;                           % Preference for smaller subspaces
            end
%             ix, size(S_{mm},1)
            condition = (ix ~= current) & abs(misa_values(ix) - misa_values(current)) < sqrt(eps);
            if all(condition)
                % Ignore values that are smallet than current if diff is
                % too small
                ix = current;
            else
                [~,tmp]=min(condition);
                ix = ix(tmp);
            end

            if ix < (ss+1)
                % Remove added subspace
                S_{mm}(ix,kk) = 1;
                S_{mm} = S_{mm}(1:ss,:);
            else
                b = [b'; b(end)];
                l = [l'; l(end)];
                e = [e'; e(end)];
            end
            % Retain only non-empty subspaces... RECONSIDER this
            this_nes = full(sum(S_{mm},2)) ~= 0;
            S_{mm} = S_{mm}(this_nes,:);
            b = b(this_nes);
            l = l(this_nes);
            e = e(this_nes);
        end
%         figure, imagesc(full(S{mm}),[0 1]), axis equal tight; title('DESIRED')
%         figure, imagesc(full(S_{mm}),[0 1]), axis equal tight; title('BEFORE subspace match: ESTIMATED')
        % Reorder components to match new S_ as close as possible to S
        shuff = match_subspaces(S{mm},S_{mm});
        S_{mm} = S_{mm}(:,shuff);
%         figure, imagesc(full(S_{mm}),[0 1]), axis equal tight; title('AFTER subspace match: ESTIMATED')
        O.update(S_,mm,b',l',e');
        % Reorder rows of W acordingly
        oldW{mm} = oldW{mm}(shuff,:);
        O.objective(O.ut.stackW({O.W{mm}(shuff,:)}));
    end
    
    O.update(S,M,old_beta,old_lambda,old_eta);
%     O.updateRElambda(REl);
    [~,shuff] = O.sub_perm_analysis(O.ut.stackW(O.W));
    O.objective(O.ut.stackW(cellfun(@(w,s) w(s,:), oldW(M), shuff(M), 'Un', 0)));
    w0 = O.ut.stackW(O.W(O.M));
    O.updatesc(sc);                                 % Set scale-control
    
else
    w0 = O.ut.stackW(O.W(O.M));
end

    function out = match_subspaces(Sdes,Sest)
        
        Sd = full(sum(Sdes,2));                     % DESired number of sources in each subspace
        [ldes,ix_] = sort(Sd(Sd ~= 0));             % Subspace sizes in ascending order
        [~,ix1] = sort(ix_);                        % Indices to sort like Sdes
        
        shuff_ = cell(1,size(Sest,1));
        for ss = 1:size(Sest,1)
            shuff_{ss} = find(Sest(ss,:));          % Source indices for ESTimated subspace
        end
        
        % if size(Sdes,1) > size(Sest,1)
        %     num_empty = size(Sdes,1) - size(Sest,1);
        % else
        %     num_empty = 0;
        % end
        % [lest,ix2] = sort([cellfun(@length,shuff) zeros(1,num_empty)]);
        [lest,ix2] = sort(cellfun(@length,shuff_)); % ESTimated subspace sizes in ascending order
        
        ld = 1; le = 1;
        missd = []; misse = [];                     % mismatched entries
        matchd = []; matche = [];                   % matched entries
        while ld <= length(ix_) && le <= length(ix2)
            if ldes(ld) == lest(le)                 % matching by subspace size
                % tmp_ix2(le) = ix2(ix(ld));
                matche = [matche le];
                matchd = [matchd ld];
                ld = ld + 1;
                le = le + 1;
            elseif ldes(ld) > lest(le)              % Mismatch: ESTimated subspace is smaller
                while (ld <= length(ix_)) && (le <= length(ix2)) && (ldes(ld) > lest(le))
                    %             if lest(le) ~= 0
                    misse = [misse le];
                    %             end
                    le = le + 1;
                end
            else                                    % Mismatch: DESired subspace is smaller
                while (ld <= length(ix_)) && (le <= length(ix2)) && (ldes(ld) < lest(le))
                    missd = [missd ld];
                    ld = ld + 1;
                end
            end
        end
        missd = [missd ld:length(ix_)];             % Append what's left
        misse = [misse le:length(ix2)];
        
        tmp_ix2 = zeros(1,max(length(ix_),length(ix2))); % largest of the 2
%         tmp = ix2(matche);         % copy matched indices
        tmp_ix2(ix_(matchd)) = ix2(matche); %tmp(ix_(matchd));         % copy matched indices
        mx = min(length(missd), length(misse));
%         tmp = ix2(misse);
        tmp_ix2(ix_(missd(1:mx))) = ix2(misse(1:mx)); %tmp(ix_(missd(1:mx))); % copy missed indices
        the_zeros = find(tmp_ix2 == 0);             % unassigned empty slots
%         tmp = ix2(misse((mx+1):end));
        tmp_ix2(the_zeros(1:length(misse((mx+1):end)))) = ...
            ix2(misse((mx+1):end)); %tmp(the_zeros(1:length(misse((mx+1):end))));                 % remaining missed go on empty slots
        tmp_ix2 = tmp_ix2(tmp_ix2 ~= 0);            % discard any still unassigned slots
        
        % list of ESTimated source indices in order that matches the
        % DESired (user provided) arrangement:
        out = [shuff_{tmp_ix2}];
        
        %ix2 = [ix2(ix(1:mx)) ix2((mx+1):end)];
        
    end

end