classdef gsd
    % Generate Simulated Data
    properties
        A                           % Matrix, mixing matrices (1 per dataset)
        C                           % Vector, number of components per dataset
        K                           % Number, number of non-empty subspaces
        M                           % Vector, indices of datasets to be used in the analysis
        N                           % Number, number of observations/measurements/samples
        S                           % Matrix, subspace assignment matrix (K x Cm)
        V                           % Vector, dimensionality of each dataset (number of sensors)
        W                           % Matrix, unmixing matrices (1 per dataset)
        Wi                          % Matrix, Wiener optimal unmixing matrices (1 per dataset)
        Y                           % Matrix, source estimates (network output)
        Ryy                         % Matrix, source covariance of each dataset
        d                           % Vector, dimensionality of each subspace
        dist                        % Struct Array, distribution parameters for each subspace
        tag                         % String, type of subspace structure
        Acond                       % Vector, condition number of A matrix
        Atype                       % String, A matrix was provided or generated
        tagA                        % CellString, square or rectangular for each dataset
        ae                          % Vector, scaling constant for noise (computed from SNR)
        SNR                         % Vector, signal-to-noise ratio
        SNRdB                       % Vector, signal-to-noise ratio in decibels (dB)
        noise                       % Matrix, additive noise (1 per dataset)
        predictors                  % Vector, indices of components used to generate ref
        regtype                     % String, 'LS' or 'Logistic'
        ref                         % Matrix, reference generated frmo columns of A
        seed                        % Struct, random number generator seed
        ut                          % Container for utility static functions
    end
    methods
        function obj = gsd(seed, N, M_Tot, K, C, S_, dist_params, V, Atype, Acond, A, SNR, predictors, regtype)
            % Set seed:
            if isempty(seed)
                rng(1981);
                obj.seed = rng;
            else
                rng(seed);
                obj.seed = rng;
            end
            
            % SOURCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.N = N;
            obj.M = 1:M_Tot;
            obj.K = K;
            obj.C = C;
            
            % Parse subspace structure information
            for mm = obj.M
                if issparse(S_{mm})
                    obj.S{mm} = S_{mm};
                else
                    ii = [];
                    jj = [];
                    for ii_ = 1:obj.K
                        jj_ = length(S_{ii_,mm});
                        if jj_ ~= 0
                            jj = [jj S_{ii_,mm}];
                            ii = [ii ii_*ones(1,jj_)];
                        end
                    end
                    obj.S{mm} = sparse(ii, jj, ones(1,obj.C(mm)), ...
                        obj.K, obj.C(mm), obj.C(mm));
                end
            end
            obj.d = full(sum([obj.S{obj.M}],2));  % Dimensionality of each subspace
            
            % Determine type of structure
            if length(obj.M) == 1
                % either ICA or ISA
                if length(unique(obj.d)) == 1
                    obj.tag = 'ICA';
                else
                    obj.tag = 'ISA';
                end
            else
                % either IVA or MISA
                if length(unique(obj.d)) == 1
                    obj.tag = 'IVA';
                else
                    obj.tag = 'MISA';
                end
            end
            
            % Generate a random sample for each SCV:
            obj.ut = utils;
            SCV = cell(1, obj.K);
            for kk = 1:obj.K
                switch(lower(dist_params(kk).name))
                    case 'mvn'
                        obj.dist(kk).name = 'Multivariate Normal';
                        
                    case 'mvl'
                        obj.dist(kk).name = 'Multivariate Laplace';
                        
                    case 'mvk'
                        obj.dist(kk).name = 'Multivariate K';
                        
                    case 'mpe'
                        obj.dist(kk).name = 'Multivariate Power Exponential';
                        
                    case 'kotz'
                        obj.dist(kk).name = 'Multivatiate Kotz';
                        
                    case 'gcp'
                        obj.dist(kk).name = 'Gaussian Copula';
                        
                end
                % Setup SCV distribution parameters:
                obj.dist(kk).mu = dist_params(kk).mu;                               % SCV mean
                if obj.d(kk) == 1
                    obj.dist(kk).r = 1;
                    obj.dist(kk).R = 1;
                else
                    options = optimset('TolX',1e-100);
                    obj.dist(kk).r =  fzero(@(a) (a)^(1/(obj.d(kk)-1)) - dist_params(kk).CORR, ...
                        [0 (1-eps)], options);                                          % max. corr.
                    obj.dist(kk).R = toeplitz(obj.dist(kk).r.^linspace(0,1,obj.d(kk))); % SCV correlation
                end
                
                switch(obj.dist(kk).name)
                    case 'Multivariate Normal'
                        SCV{kk} = mvnrnd(obj.dist(kk).mu',obj.dist(kk).R,obj.N);
                        
                    case 'Multivariate Laplace'
                        nes = obj.d~=0;
                        eta = 1;
                        beta = .5;
                        nu = exp(log(2*eta + obj.d(kk) - 2) - log(2*beta));
                        desired_stddev = pi/sqrt(3);
                        obj.dist(kk).lambda = exp(beta * (gammaln(nu + 1./beta) - log(obj.d(kk)) - gammaln(nu) - 2*log(desired_stddev)));
                        
                        %obj.dist(kk).lambda = .5;
                        
                        % Adjusted to be stddev = 1:
                        SCV{kk} = ...sqrt((1/(obj.d(kk)+1)) * det(obj.dist(kk).R)^(1/obj.d(kk))) * ...
                            obj.ut.mymvlap(obj.dist(kk).mu,obj.dist(kk).R,obj.N,1/obj.dist(kk).lambda);
                        
                    case 'Multivariate K'
                        if length(dist_params.a) == 1
                            obj.dist(kk).a = dist_params.a*ones(1, obj.K);
                        else
                            obj.dist(kk).a = dist_params.a;
                        end
                        if length(dist_params.b) == 1
                            obj.dist(kk).b = dist_params.b*ones(1, obj.K);
                        else
                            obj.dist(kk).b = dist_params.b;
                        end
                        
                        % Adjusted to be stddev = 1:
                        SCV{kk} = ...sqrt((1/(obj.d(kk)+1)) * det(obj.dist(kk).R)^(1/obj.d(kk))) * ...
                            obj.ut.mymvk(obj.dist(kk).mu,obj.dist(kk).R,obj.N,obj.dist(kk).a,obj.dist(kk).b);
                        
                    case 'Multivariate Power Exponential'
                        
                    case 'Multivatiate Kotz'
                        
                    case 'Gaussian Copula'
                        U = copularnd('Gaussian',obj.dist(kk).R,obj.N);
                        out = zeros(size(U));
                        alpha = sqrt(size(U,2) + 1);
                        for cc = 1:size(U,2)
                            % Univariate paramters
                            out(:,cc) = obj.ut.myicdf('laplace', U(:,cc), obj.dist(kk).mu(cc), alpha);
                        end
                        SCV{kk} = out;
                        
                end
            end
            obj.dist = obj.dist';
            
            % Assign SCVs to datasets
            obj.Y = cell(1,max(obj.M));

            for kk = 1:obj.K
                ss = 0;
                set_ = 1:size(SCV{kk},2);
                done = false(1,max(obj.M));
                selected = cell(1,max(obj.M));
                
                while ~all(done)
                    skip = max(obj.M);
                    for mm = obj.M
                        %obj.Y{obj.M(mm)} = [];
                        %obj.Y{obj.M(mm)} = [obj.Y{obj.M(mm)}; SCV{kk}(:,ss+(1:ss_))'];
                        if ~done(mm)
                            ss_ = sum(logical(obj.S{mm}(kk,:)));
                            how_many = ss_ - length(selected{mm});
                            
                            [selected_, done_] = select_from_set(obj, set_, how_many, skip);
                            done(mm) = done_;
                            selected{mm} = [selected{mm}, selected_];
                            set_ = setdiff(set_, selected{mm});
                            skip = skip - 1;
                        end
                    end
                end
                
                for mm = obj.M
                    ss_ = sum(logical(obj.S{mm}(kk,:)));
                    numY = size(obj.Y{obj.M(mm)},1);
                    obj.Y{obj.M(mm)} = [obj.Y{obj.M(mm)}; SCV{kk}(:,selected{mm})'];
                    if kk == 1
                        obj.Ryy{mm} = zeros(obj.C(mm));
                    end
                    if strcmpi(obj.dist(kk).name, 'Multivariate K') || ...
                        strcmpi(obj.dist(kk).name, 'Gaussian Copula')
                        obj.Ryy{mm}(numY+(1:ss_),numY+(1:ss_)) = cov(obj.Y{obj.M(mm)}(numY+(1:ss_),:)');
                        % obj.Ryy{mm}(selected{mm},selected{mm}) = cov(obj.Y{obj.M(mm)}(selected{mm},:)');
                    else
                        obj.Ryy{mm}(numY+(1:ss_),numY+(1:ss_)) = (obj.d(kk)+1) * ...
                            obj.dist(kk).R(ss+(1:ss_),ss+(1:ss_));
                        % obj.Ryy{mm}(selected{mm},selected{mm}) = (obj.d(kk)+1) * ...
                        %     obj.dist(kk).R(selected{mm},selected{mm});
                    end
                    ss = ss + ss_;
                end
            end
            
            % MIXTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Quadratic? Non-Lin?
            
            obj.Atype = Atype;
            switch (obj.Atype)
                case 'Provided'
                    obj.A = A;
                    for mm = obj.M
                        Acond = cond(A{mm});
                        obj.V(mm) = size(obj.A{mm},1);
                        obj.Acond(mm) = Acond*ones(size(mm));
                        obj.W{mm} = pinv(A{mm});
                    end
                case 'Generated'
                    if isempty(Acond)
                        obj.Acond(obj.M) = 3;
                    else
                        %assume double
                        if length(Acond) == 1
                            obj.Acond(obj.M) = Acond*ones(size(obj.M));
                        else
                            obj.Acond = Acond;
                        end
                    end
                    
                    if isempty(V) && isempty(A)
                        obj.V = obj.C; % default to square mixing
                    else
                        % assume double
                        if length(V) == 1
                            obj.V(obj.M) = V*ones(size(obj.M));
                        else
                            obj.V = V;
                        end
                    end
                    
                    obj.A = cell(1,max(obj.M));
                    for mm = obj.M
                        [u,s,v] = svd(randn(obj.V(mm),obj.C(mm)), 'econ');
                        if obj.Acond(mm) == 1
                            obj.A{mm} = u*v';
                            obj.W{mm} = v*u';
                        else
                            obj.A{mm} = bsxfun(@times,u,diag(s)'+((max(diag(s))-obj.Acond(mm)*min(diag(s)))/(obj.Acond(mm)-1)))*v';
                            obj.W{mm} = bsxfun(@times,v,1./(diag(s)'+((max(diag(s))-obj.Acond(mm)*min(diag(s)))/(obj.Acond(mm)-1))))*u';
                        end
                    end
                    
                    for mm = obj.M
                        if obj.V(mm) == obj.C(mm)
                            obj.tagA{mm} = 'square';
                        else
                            obj.tagA{mm} = 'rectangular';
                        end
                    end
                    
            end
            
            % NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is white noise. Otherwise should use mvnrnd, define Rnn
            % and replace obj.V by trace(Rnn).
            if isempty(SNR)
                obj.SNR = (999+1)/1; % default to SNRdB = 30dB
            else
                % assume double
                if length(SNR) == 1
                    obj.SNR(obj.M) = SNR*ones(size(obj.M));
                else
                    obj.SNR = SNR;
                end
            end
            
            for mm = obj.M
                %obj.ae(mm) = sqrt((trace(obj.A{mm}*obj.Ryy{mm}*obj.A{mm}')/obj.V(mm))/(obj.SNR(mm)-1));
                obj.ae(mm) = sqrt((sum(sum(obj.A{mm}.*(obj.A{mm}*obj.Ryy{mm}),2))/obj.V(mm))/(obj.SNR(mm)-1));
                obj.noise{mm} = sqrt(obj.ae(mm))*randn(obj.V(mm),obj.N);
            end
            obj.SNRdB = 10*log10(obj.SNR);
            
            % WIENER OPTIMAL UNMIXING (NOISY) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is for white noise. Otherwise should replace 1/V by
            % Rnn^-1
            for mm = obj.M
                if obj.ae(mm) ~= 0
                    obj.Wi{mm} = ((1/obj.ae(mm))*(obj.A{mm}'*obj.A{mm}) + inv(obj.Ryy{mm}))\(obj.A{mm}'*(1/obj.ae(mm)));
                else
                    obj.Wi{mm} = obj.W{mm};
                end
            end
            
            % REFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if iscell(predictors)
                obj.predictors = predictors;
            else
                % assume double
                for mm = obj.M
                    obj.predictors{mm} = predictors;%{[1:2], [], [1:2]};
                end
            end
            
            obj.regtype = regtype;
            if isempty(obj.regtype)
                if isempty([obj.predictors{:}])
                    obj.regtype = '';
                else
                    obj.regtype = 'LS';
                end
            end
            
            
            switch obj.regtype
                case 'LS'
                    beta = [5 5]';
                    for mm = obj.M
                        if isempty(obj.predictors{mm})
                            obj.ref{mm} = [];
                        else
                            obj.ref{mm} = obj.A{mm}(:,obj.predictors{mm})*beta + ...
                                .01*randn(size(obj.A{mm},1),1);
                        end
                    end
                case 'Logistic'
                    beta = [5 10]';
                    for mm = obj.M
                        if isempty(obj.predictors{mm})
                            obj.ref{mm} = [];
                        else
                            pr = 1./(1+exp(-obj.A{mm}(:,obj.predictors{mm})*beta));
                            obj.ref{mm} = random('Binomial',1,pr);
                        end
                    end
                otherwise
                    obj.ref = {};
            end
            
        end
        X = genX(O);                    % Generate the data X
        saveX(O,path,prefix,suffix);    % Save object to .mat file
        saveme(O,path,prefix,suffix);   % Save object to .mat file
        [selected, done] = select_from_set(O, set_, how_many, skip);
    end
    
end