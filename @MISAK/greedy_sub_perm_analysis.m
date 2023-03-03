function [shuff] = greedy_sub_perm_analysis(O)

% figure,imagesc(O.S{1},max(max(abs(O.S{1}))).*[-1 1]);colorbar();
% figure,imagesc(O.S{end},max(max(abs(O.S{end}))).*[-1 1]);colorbar();

S_ = O.S;
col_ind = cell(1, O.K);
for ss = 1:O.K % iterate through subspaces
    col_ind{ss} = cellfun(@(s) find(s(ss,:)), S_, 'un', 0);
end

for tt = randperm(O.K) % iterate through subspaces
    % find the non-empty modality, get the first non-empty column
    non_empty_ind = cellfun( @(c) ~isempty(c), col_ind{tt} );
    mm = find(non_empty_ind);
    mm = mm(1);
    cc = col_ind{tt}{O.M(mm)}(1);
    current = find(S_{O.M(mm)}(:,cc)); % current subspace for component cc

    % get components from current subspace per modality: cell array index
    % for each modality
    % remove components from current subspace
    kk = cellfun(@(s) find(s(current,:)), S_, 'un', 0);

    % make index to 0
    for mm = O.M
        S_{mm}(:,kk{mm}) = zeros(size(S_{mm},1),length(kk{mm}));   % empty column cc
    end

    misa_values = [];                         % container for cost values

    % iterate through 1 to find(O.nes)
    % replace mm with for-loop over modalities
    ind_nes_subspace = find(O.nes);
    for ss = ind_nes_subspace
        ind_nes_ss1 = ind_nes_subspace(1);
        if ss ~= ind_nes_ss1
            for mm = O.M
                S_{mm}(ss_prev,kk{mm}) = 0;   % Remove previous assignment
            end
        end
        for mm = O.M
            S_{mm}(ss,kk{mm}) = 1;            % Assign component cc to subspace ss
        end
        % update on M
        O.update(S_,O.M,O.beta(ind_nes_ss1),[],O.eta(ind_nes_ss1));
        w0 = O.ut.stackW(O.W);
        misa_values(ss) = O.objective(w0);    % Compute cost
        ss_prev = ss;
    end

    for mm = O.M
        S_{mm}(ss_prev,kk{mm}) = 0;           % Remove previous assignment
    end

    [~,ix] = min(misa_values);                % best subspace for component cc
    %             imagesc(full(S_{mm}))
    drawnow()
    ix = find(misa_values == misa_values(ix));

    condition = (ix ~= current) & abs(misa_values(ix) - misa_values(current)) < sqrt(eps);
    if all(condition)
        % Ignore values that are smaller than current if diff is
        % too small
        ix = current;
    else
        [~,tmp]=min(condition);
        ix = ix(tmp);
    end

    for mm = O.M
        S_{mm}(ix,kk{mm}) = 1;
    end

    % Retain only non-empty subspaces... RECONSIDER this
    this_nes = full(sum([S_{:}],2)) ~= 0;
    S_ = cellfun(@(s) s(this_nes,:), S_, 'un', 0);
    b = O.beta(this_nes);
    e = O.eta(this_nes);

end

O.update(S_,O.M,b(1),[],e(1));

shuff = cell(1,max(O.M));
for mm = O.M
    shuff{mm} = [];
    for ss = 1:size(S_{O.M(1)},1)                   % Loop through suspaces
        shuff{mm} = [shuff{mm} find(S_{mm}(ss,:))];
    end
end

% Apply shuffle to original structure
O.updatesc(true);
O.update(O.S,O.M,O.beta(1),[],O.eta(1));
O.objective(O.ut.stackW(cellfun(@(w,s) w(s,:), O.W(O.M), shuff(O.M), 'Un', 0)));

end