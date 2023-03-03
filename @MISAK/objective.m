function [J, gJ] = objective(O,w)
 
O.W(O.M) = O.ut.unstackW(w,O.M,O.C,O.V);
O.Y(O.M) = cellfun(@mtimes, O.W(O.M), O.X(O.M), 'Un', 0);
if O.sc == false
    if nargout > 1
        [J, gJ] = objective_(O);
    else
        J = objective_(O);
    end
else
    if nargout > 1
        [J, gJ] = objective_sc_(O);
    else
        J = objective_sc_(O);
    end
end

if nargout > 1
    gJ = O.ut.stackW(gJ(O.M));
end

end