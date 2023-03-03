function setREref(O,REr,RErtype,RErlambda,rC)

if ~isempty(REr)
    O.REref = REr;
    O.ner = false(1,max(O.M));
    for mm = O.M
        if ~isempty(O.REref{mm})
            O.ner(mm) = true;
        end
    end
    if ~isempty(RErtype)
        if ischar(RErtype)
            switch lower(RErtype)
                case 'linearreg'
                    O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                    O.REreftype = cell(1,max(O.M));
                    O.REreftype(O.M) = cellstr(repmat('linearreg',length(O.M),1)); % linear regression with intercept
                case 'linearreg_nointercept'
                    O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                    O.REreftype = cell(1,max(O.M));
                    O.REreftype(O.M) = cellstr(repmat('linearreg_nointercept',length(O.M),1)); % linear regression without intercept
                case 'logisticls'
                    O.REref(O.M) = cellfun(@(x) log((x+.5)./(1-x+.5)),O.REref,'Un',0);
                    O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                    O.REreftype = cell(1,max(O.M));
                    O.REreftype(O.M) = cellstr(repmat('logisticls',length(O.M),1)); % logistic regression by logit least squares
                case 'logisticirls'
                    O.REref(O.M) = cellfun(@(x) log((x+.5)./(1-x+.5)),O.REref,'Un',0);
                    O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                    O.REreftype = cell(1,max(O.M));
                    O.REreftype(O.M) = cellstr(repmat('logisticirls',length(O.M),1)); % logistic regression by IRLS
                otherwise
                    warning('Unknown reference type (REreftype): %s.\nSwitching to default: linearreg', RErtype)
                    O.REreftype = cell(1,max(O.M));
                    O.REreftype(O.M) = cellstr(repmat('linearreg',length(O.M),1)); % linear regression with intercept
            end
        elseif iscell(RErtype)
            if length(RErtype) == 1
                switch lower(RErtype{1})
                    case 'linearreg'
                        O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                        O.REreftype = cell(1,max(O.M));
                        O.REreftype(O.M) = cellstr(repmat('linearreg',length(O.M),1)); % linear regression with intercept
                    case 'linearreg_nointercept'
                        O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                        O.REreftype = cell(1,max(O.M));
                        O.REreftype(O.M) = cellstr(repmat('linearreg_nointercept',length(O.M),1)); % linear regression without intercept
                    case 'logisticls'
                        O.REref(O.M) = cellfun(@(x) log((x+.5)./(1-x+.5)),O.REref,'Un',0);
                        O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                        O.REreftype = cell(1,max(O.M));
                        O.REreftype(O.M) = cellstr(repmat('logisticls',length(O.M),1)); % logistic regression by logit least squares
                    case 'logisticirls'
                        O.REref(O.M) = cellfun(@(x) log((x+.5)./(1-x+.5)),O.REref,'Un',0);
                        O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
                        O.REreftype = cell(1,max(O.M));
                        O.REreftype(O.M) = cellstr(repmat('logisticirls',length(O.M),1)); % logistic regression by IRLS
                    otherwise
                        warning('Unknown reference type (REreftype): %s.\nSwitching to default: linearreg', RErtype)
                        O.REreftype = cell(1,max(O.M));
                        O.REreftype(O.M) = cellstr(repmat('linearreg',length(O.M),1)); % linear regression with intercept
                end
            else
                O.REreftype = RErtype;
            end
        else
            error('REreftype must be either cell or char.');
        end
    else
        O.RErefnorm = cellfun(@(x) mean(sum(x.^2,2)),O.REref,'Un',0);
        O.REreftype = cell(1,max(O.M));
        O.REreftype(O.M) = cellstr(repmat('linearreg',length(O.M),1)); % linear regression with intercept
    end
    if ~isempty(RErlambda)
        if isfloat(RErlambda)
            if length(RErlambda) == 1
                O.REreflambda = cell(1,max(O.M));
                O.REreflambda(O.ner) = {RErlambda}; % linear regression with intercept
            else
                O.REreflambda = mat2cell(RErlambda,1,length(RErlambda));
            end
        elseif iscell(RErlambda)
            if length(RErlambda) == 1
                O.REreflambda = cell(1,max(O.M));
                O.REreflambda(O.ner) = RErlambda; % linear regression with intercept
            else
                O.REreflambda = RErlambda;
            end
        else
            error('REreflambda must be either cell or double.');
        end
        
    else
        O.REreflambda = cell(1,max(O.M));
        O.REreflambda(O.M) = {.5}; % proportion of unexplained variance
    end
    if ~isempty(rC)
        if isfloat(rC)
            if length(rC) == 1
                O.rC = cell(1,max(O.M));
                O.rC(O.ner) = {rC};
            else
                O.rC = mat2cell(rC,1,length(rC));
            end
        elseif iscell(rC)
            if length(rC) == 1
                O.rC = cell(1,max(O.M));
                O.rC(O.ner) = rC;
            else
                O.rC = rC;
            end
        else
            error('rC must be either cell or double.')
        end
    else
        O.rC = cell(1,max(O.M));
        O.rC(O.M) = mat2cell(O.C, 1, length(O.C));
    end
    %     if ~isempty(rC)
    %         if isfloat(rC)
    %             if length(rC) == 1
    %             else
    %             end
    %         elseif iscell(rC)
    %             if length(rC) == 1
    %             else
    %             end
    %         end
    %     else
    %         error('rC must be either cell or double.')
    %     end
else
    O.REref = {};
    O.RErefnorm = {};
    O.REreftype = {''};
    O.REreflambda = {};
    O.rC = {};
end

end