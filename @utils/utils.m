classdef utils
    properties
    end
    methods
        function obj = utils(s)
        end
        [ISI, WAr] = MISI(O,W,A,S)
        [MD, WAr] = MMD(O,W,A,S)
        [MSE, WAr] = MMSE(O,W,Y,Yh)
        MD = myMD(O,WAr)
        MSE = myMSE(O,WAr)
        [out, R_, S] = mymvlap(O,mu,R,N,b)
    end
    %     methods (Access = private)
    %         output = myFunc(obj,arg1,arg2)
    %     end
    methods (Static)
        w = stackW(W)
        W = unstackW(w,M,C,V)
        w = stackMuCov(mu_, cov_)
        [mu_, cov_] = unstackMuCov(w,M,V)
        ISI = myISI(WAr)
        [assignment,cost] = munkres(WAr)
        %         Vt = myorth(Y)
        [out, R_, S] = mymvk(mu,R,N,a,b)
        optprob = getop(w0,varargin)
        X = myicdf(name,U,varargin)
        [Ywht,whtM,dewhtM,Ybr,brM,debrM] = myPCA(Y)
        [whtM, H] = doMMGPCA(X, comps, rec_type)
        [wout,fval,exitflag,output] = run_MISA(misa_obj,W0)
    end
end % End of classdef

% function myUtilityFcn
% end