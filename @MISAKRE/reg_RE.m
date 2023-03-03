function [c,DC] = reg_RE(O,w)

% check if need to update W
if sum(O.ut.stackW(O.W) ~= w) > 0
    O.W(O.M) = O.ut.unstackW(w,O.M,O.C,O.V);
    O.Y(O.M) = cellfun(@mtimes, O.W(O.M), O.X(O.M), 'Un', 0);
end

if nargout > 1
    [c, DC] = O.objective(w);
    [c_, DC_] = O.opt_RE(w);
    c = c + (1e-2/O.RElambda)*(c_ + O.RElambda);%(c_ + O.RElambda);%
    DC = DC + (1e-2/O.RElambda)*(DC_ + O.RElambda);%(DC_ + O.RElambda);%
else
    c = O.objective(w) + (1/O.RElambda)*(O.opt_RE(w) + O.RElambda);
end

end