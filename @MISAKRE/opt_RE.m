function [c,DC] = opt_RE(O,w)

if sum(O.ut.stackW(O.W(O.M)) ~= w) > 0
    O.W(O.M) = O.ut.unstackW(w,O.M,O.C,O.V);
    O.Y(O.M) = cellfun(@mtimes, O.W(O.M), O.X(O.M), 'Un', 0);
end

if nargout > 1
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
                        DC_(ii) = DC_(ii) + O.ut.stackW(DC(2,O.ner & sel));
                    end
                elseif ~isempty(O.REref)
                    auxmsk = O.ner & sel;
                    for mm = O.M
                        if auxmsk(mm)
                            ii = [ii; sum(O.C(O.M(1:(mm-1))).*O.V(O.M(1:(mm-1)))) + (1:(O.C(mm)*O.V(mm)))'];
                        end
                    end
                    DC_ = zeros(size(O.ut.stackW(DC(1,O.M))));
                    DC_(ii) = O.ut.stackW(DC(2,O.ner & sel));
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
                        DC_(ii) = DC_(ii) + O.ut.stackW(DC(2,O.ner & sel));
                    end
                elseif ~isempty(O.REref)
                    auxmsk = O.ner & sel;
                    for mm = O.M
                        if auxmsk(mm)
                            ii = [ii; sum(O.C(O.M(1:(mm-1))).*O.V(O.M(1:(mm-1)))) + (1:(O.C(mm)*O.V(mm)))'];
                        end
                    end
                    DC_ = zeros(size(O.ut.stackW(DC(1,O.M))));
                    DC_(ii) = O.ut.stackW(DC(2,O.ner & sel));
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
        end
    else
        DC = [];
    end
else
    c = O.RE_();
end

if size(c,1) > 1
    c = sum([([c{1,:}] - O.RElambda) ([c{2,:}] - [O.REreflambda])]);
else
    c = sum([c{1,:}]) - O.RElambda;
end

end

% if nargout > 1
%     [c, DC] = O.RE_();
% else
%     c = O.RE_();
% end
% 
% c = sum([c{:}] + 1e-20);
% 
% if nargout > 1
%     switch O.REtype_
%         case true                       % 'NMSE' (log of)
%             for mm = O.M
%                 DC{mm} = (2/(c*O.Xnorm(mm)*log(10)))*DC{mm};
%             end
%             DC = O.stackW(DC(O.M));
%         case false                      % 'MSE'
%             DC = O.stackW(DC(O.M));
%             DC = (2/(c*log(10)))*DC;
%     end
% end
% 
% c = (log(c) - log(O.N))/log(10) - O.RElambda;
