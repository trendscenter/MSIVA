function [c,ceq,DC,DCeq] = con_RE(O,w)

c = [];
ceq = [];

if sum(O.ut.stackW(O.W(O.M)) ~= w) > 0
    O.W(O.M) = O.ut.unstackW(w,O.M,O.C,O.V);
    O.Y(O.M) = cellfun(@mtimes, O.W(O.M), O.X(O.M), 'Un', 0);
end

if nargout > 2
    [c, DC] = O.RE_();
    if ~isempty(DC)
        sel = false(1, max(O.M));
        sel(O.M) = true; % selected datasets
        if length(O.M) == 1
            if size(c,1) > 1
                ii = [];
                DC_ = [];
                if ~isempty(O.REtype_)
                    DC_ = O.ut.stackW(DC(1,O.M));
                    if ~isempty(O.REref)
                        auxmsk = O.ner & sel;
                        for mm = O.M
                            if auxmsk(mm)
                                ii = [ii; sum(O.C(O.M(1:(mm-1))).*O.V(O.M(1:(mm-1)))) + (1:(O.C(mm)*O.V(mm)))'];
                            end
                        end
                        DC_(ii,2) = O.ut.stackW(DC(2,O.ner & sel));
                    end
                elseif ~isempty(O.REref)
                    auxmsk = O.ner & sel;
                    for mm = O.M
                        if auxmsk(mm)
                            ii = [ii; sum(O.C(O.M(1:(mm-1))).*O.V(O.M(1:(mm-1)))) + (1:(O.C(mm)*O.V(mm)))'];
                        end
                    end
                    DC_ = zeros(size(O.ut.stackW(DC(1,O.M))));
                    DC_(ii,1) = O.ut.stackW(DC(2,O.ner & sel));
                end
                DC = DC_;
                clear DC_
            else
                if ~isempty(O.REtype_)
                    DC = O.ut.stackW(DC(1,O.M));
                else
                    DC = [];
                end
            end
        else
            if size(c,1) > 1
                ii = [];
                jj = [];
                if ~isempty(O.REtype_)
                    tot = sum(O.C(O.M).*O.V(O.M));
                    ii = (1:tot)';
                    jj = O.ut.stackW(cellfun(@(mm,x) bsxfun(@times,mm,ones(size(x))), ...
                        num2cell(1:length(O.M)), O.W(O.M), 'Un', 0));
                    if ~isempty(O.REref)
                        nerTotal = sum(O.ner(O.M));
                        auxmsk = O.ner & sel;
%                         auxind = find(auxmsk);
                        for mm = O.M
                            if auxmsk(mm)
                                ii = [ii; sum(O.C(O.M(1:(mm-1))).*O.V(O.M(1:(mm-1)))) + (1:(O.C(mm)*O.V(mm)))'];
                            end
                        end
%                         newind = [];
%                         for ff = 1:length(auxind)
%                             newind = [newind find(O.M == auxind(ff))];
%                         end
                        jj = [jj; O.ut.stackW(cellfun(@(mm,x) bsxfun(@times,mm,ones(size(x))), ...
                            num2cell(length(O.M)+(1:nerTotal)), O.W(auxmsk), 'Un', 0))];
                    end
                    DC = sparse(ii, jj, [O.ut.stackW(DC(1,O.M))' O.ut.stackW(DC(2,O.ner & sel))'], ...
                        tot, length(O.M)+nerTotal, tot+sum(O.C(O.ner & sel).*O.V(O.ner & sel)));
                elseif ~isempty(O.REref)
                    nerTotal = sum(O.ner(O.M));
                    auxmsk = O.ner & sel;
                    for mm = O.M
                        if auxmsk(mm)
                            ii = [ii; sum(O.C(O.M(1:(mm-1))).*O.V(O.M(1:(mm-1)))) + (1:(O.C(mm)*O.V(mm)))'];
                        end
                    end
                    jj = [jj; O.ut.stackW(cellfun(@(mm,x) bsxfun(@times,mm,ones(size(x))), ...
                        num2cell(1:nerTotal), O.W(auxmsk), 'Un', 0))];
                    DC = sparse(ii, jj, O.ut.stackW(DC(2,O.ner & sel)), ...
                        tot, nerTotal, sum(O.C(O.ner & sel).*O.V(O.ner & sel)));
                else
                    DC = [];
                end
            else
                if ~isempty(O.REtype_)
                    tot = sum(O.C(O.M).*O.V(O.M));
                    ii = 1:tot;
                    jj = O.ut.stackW(cellfun(@(mm,x) bsxfun(@times,mm,ones(size(x))), ...
                        num2cell(1:length(O.M)), O.W(O.M), 'Un', 0));
                    DC = sparse(ii, jj, O.ut.stackW(DC(1,O.M)), tot, length(O.M), tot);
                else
                    DC = [];
                end
            end
        end
        DCeq = [];
    else
        DC = [];
        DCeq = [];
    end
else
    c = O.RE_();
end

if ~isempty(c)
    if size(c,1) > 1
        c = [([c{1,:}] - O.RElambda) ([c{2,:}] - [O.REreflambda{:}])];
    else
        c = [c{1,:}] - O.RElambda;
    end
else
    c = [];
end

end