classdef MISAKRE < MISAK
    properties
        Xnorm
        REtype      % NMSE/MSE
        REapproach  % W' / pinv(W)
        RElambda
        REref
        RErefnorm
        REreftype
        REreflambda
        rC          % list of component indexes used to estimate reference
        ner         % non-empty reference
    end
    properties (Access = private)
        REtype_
        REapproach_
    end
    methods
        function obj = MISAKRE(w0, M, S, X, beta, eta, lambda, gradtype, sc, preX, ...
                REtype, REapproach, RElambda, REref, REreftype, REreflambda, rC)
            obj@MISAK(w0, M, S, X, beta, eta, lambda, gradtype, sc, preX);
            obj.Xnorm = cell2mat(cellfun(@(x) mean(sum(x.^2)),obj.X,'Un',0));
            setREtype(obj,REtype)
            setREapproach(obj,REapproach)
            updateRElambda(obj,RElambda)
            setREref(obj,REref,REreftype,REreflambda,rC)
        end
%         [J, gJ] = objective_(O)
        [c,DC] = RE_(O)
        [c,ceq,DC,DCeq] = con_RE(O,w)
        [c,DC] = opt_RE(O,w)
        [c,DC] = reg_RE(O,w)
        setREtype(O,REt)
        setREapproach(O,REa)
        updateRElambda(O,RElambda)
        setREref(obj,REref,REreftype,REreflambda,rC)
        updateREreflambda(O,REreflambda)
    end
    methods (Access = private)
        Vt = orthY(O)
    end
%     methods (Static)
%         w = stackW(W)
%         W = unstackW(w,M,C,V)
%     end
end % End of classdef

% function myUtilityFcn
% end