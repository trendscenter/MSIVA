classdef MISAK < handle
    % MISA using Kotz distribution and intrinsic source scaling
    properties
        C                           % Vector, number of components per dataset
        K                           % Number, number of non-empty subspaces
        M                           % Vector, indices of datasets to be used in the analysis
        N                           % Number, number of observations/measurements/samples
        S                           % Matrix, subspace assignment matrix (K x Cm)
        V                           % Vector, dimensionality of each dataset (number of sensors)
        W                           % Matrix, unmixing matrices (1 per dataset)
        X                           % Matrix, all datasets
        Y                           % Matrix, source estimates (network output)
        beta                        % Vector, Kotz beta parameter for each subspace
        d                           % Vector, dimensionality of each subspace
        d_k                         % Vector, dimensionality of each subspace per dataset
        eta                         % Vector, Kotz eta parameter for each subspace
        gradtype                    % String, type of gradient used: regular, relative, or natural
        sc                          % Logical, scale-control option: true or false
        lambda                      % Vector, Kotz lambda parameter for each subspace
        nes                         % Vector, non-empty subspace indexes
        nu                          % Vector, Kotz nu parameter for each subspace
        a                           % Vector, scaling factor to convers inverse cov. into inverse dispersion mtx.
        preX                        % Logical, preprocessing option: true, false, or empty
    end
    properties (Access = protected)
        ut                          % Container for utility static functions
    end
    methods
        function obj = MISAK(w0, M, S, X, beta, eta, lambda, gradtype, sc, preX)
            
            %%% Test nu for all subspaces!!
            Sfull = full(sum([S{M}],2));
            Smask = Sfull~=0;
            if any(beta(Smask) <= 0), error('All beta parameters should be positive.'); end
            if ~isempty(lambda)
                if any(lambda(Smask) <= 0), error('All lambda parameters should be positive.'); end
            end
            if any(eta(Smask) <= ((2-Sfull(Smask))./2)), error('All eta parameters should be lagerer than (2-d)/2.'); end
            nu = (2*eta + Sfull - 2)./(2*beta);
            if any(nu(Smask) <= 0), error('All nu parameter derived from eta and d should be positive.'); end
            
            obj.M = M;                      % Indices of datasets to be used in the analysis
            obj.S = S;                      % Subspace assignment matrix (K x Cm)
            obj.X = X;                      % All datasets
            obj.beta = beta;                % Kotz beta parameter
            obj.eta = eta;                  % Kotz eta parameter
            obj.lambda = lambda;            % Kotz lambda parameter
            obj.nu = nu;                    % Kotz nu parameter
            
            if ~isempty(preX) && preX == true
                % Remove mean from all datasets:
                obj.preX.mean = cellfun(@(x) mean(x,2), obj.X, 'Un', 0);
                obj.X = cellfun(@(x,mx) bsxfun(@minus,x,mx), obj.X, obj.preX.mean, 'Un', 0);
                
                % Standardize all datasets:
                obj.preX.std = cellfun(@(x) std(x,[],2), obj.X, 'Un', 0);
                obj.X = cellfun(@(sx,x) bsxfun(@times,1./sx,x), obj.preX.std, obj.X, 'Un', 0);
                
            else obj.preX = false;
            end
            
            obj.C = cellfun(@(s) size(s,2), obj.S); % Number of components per dataset
            obj.V = cellfun(@(x) size(x,1), obj.X); % Dimensionality of each dataset (number of sensors)
            
            obj.ut = utils;
            
            % Unmixing matrices (1 per dataset):
            obj.W = cell(1,max(obj.M));
            obj.W(obj.M) = obj.unstackW(w0);%, obj.M, obj.C(obj.M), obj.V(obj.M));
            
            % Compute source estimates (network output):
            obj.Y = cell(1,max(obj.M));
            obj.Y(obj.M) = cellfun(@mtimes, obj.W(obj.M), obj.X(obj.M), 'Un', 0);
            
            obj.N = size(obj.X{obj.M(1)},2);% Number of samples
            
            obj.d = full(sum([obj.S{obj.M}],2));  % Dimensionality of each subspace
            obj.nes = obj.d~=0;             % Non-empty subspace indexes
            obj.d_k = cellfun(@(s) full(sum(s,2)), obj.S,'Un',0);
            obj.auto_tune('lambda', lambda);
            
            obj.a = (obj.lambda.^(-1./(obj.beta)) .* gamma(obj.nu + 1./obj.beta)) ./ ...
                (obj.d .* gamma(obj.nu));
            
            obj.K = sum(obj.nes);           % Number of non-empty subspaces
            if ~isempty(sc)
                updatesc(obj,sc);           % Determines if scale-control is used: true or false
            else
                updatesc(obj,true);         % Determines if scale-control is used: true or false
            end
            updategradtype(obj,gradtype);   % Determines the type of gradient used: regular, relative, natural
            
        end
        [J, gJ] = objective_(O)         % Compute objective function value at current W
        [J, gJ] = objective_sc_(O)      % Compute scale-control objective function value at current W
        [J, gJ] = objective(O,w)        % Compute objective function value at provided W
        W = unstackW(O,w)
        w = stackW(O,W)
        update(O,S,M,b,l,e)             % Subset the data
        Sold = IVAfy(O)                        % Turn problem into Unidimensional-type
        updateCS(O,S)
        updategradtype(O,gradtype)      % Set gradient type
        combinatorial_optim(O,myM)      % Solve combinatorial optimization
        w0 = greedysearch(O,myM)
        [w0, shuff] = sub_perm_analysis(O, w0)
        [shuff] = greedy_sub_perm_analysis(O)
        mISI = MISI(O,A)                % Compute joint ISI of current W. A is provided.
        mmd = MMD(O,A)
        mmse = MMSE(O,Y)
    end
    methods (Access = private)
%         output = myFunc(obj,arg1,arg2)
    end
    methods (Static)
%         w = stackW(W)
%         W = unstackW(w,M,C,V)
%         Vt = myorth(Y)
    end
end % End of classdef

% function myUtilityFcn
% end