function auto_tune(O, param, value)

if strcmpi(param,'eta')
    eta = value;
    if isempty(eta)
        eta = ones(length(O.nes),1);%ones(size(O.S{O.M(1)},1),1);
    end
    if length(eta) ~= length(O.nes)%size(O.S{O.M(1)},1)
        eta = eta(1).*ones(length(O.nes),1);%ones(size(O.S{O.M(1)},1),1);
    end
    min_val = (2-sum([O.S{O.M}],2))./2;
    if any(eta(O.nes) <= min_val(O.nes)), error('All eta parameters should be lagerer than (2-d)/2.'); end
    O.eta = eta;                  % Kotz eta parameter
    
    return
end

if strcmpi(param,'beta')
    beta = value;
    if isempty(beta)
        % beta = gammaln(pi)*3/4; % 0.617263173999417;
        % Get kurtosis-secific beta
        % 3 * (Formula 26) in this paper (3* is to scale to typical excess-kurtosis scale):
        % K. Zografos (2008), On Mardia�s and Song�s measures of kurtosis in elliptical distributions
        beta = zeros(length(O.nes),1);%zeros(size(O.S{O.M(1)},1),1);
        desired_kurtosis = 1.2; %1.591131166811146; % value that matches excess-kurtosis of logistic: ex-k = 1.2
        options = optimset('TolX',sqrt(eps));
        for kk = find(O.nes')
            %f = @(bt) 3*((((O.d(kk)^2)./(O.d(kk)*(O.d(kk)+2))) * ((gamma((2*O.eta(kk)+O.d(kk)+2)/(2*bt))*gamma((2*O.eta(kk)+O.d(kk)-2)/(2*bt)))/(gamma((2*O.eta(kk)+O.d(kk))/(2*bt))^2))) - 1) - desired_kurtosis;
            t1 = 2*log(O.d(kk))-log(O.d(kk))-log(O.d(kk)+2);
            g = @(bt) 3*(exp(t1 + (gammaln((2*O.eta(kk)+O.d(kk)+2)/(2*bt)) + gammaln((2*O.eta(kk)+O.d(kk)-2)/(2*bt)) - 2*gammaln((2*O.eta(kk)+O.d(kk))/(2*bt)))) - 1) - desired_kurtosis;
            beta(kk) = fzero(g,[7e4*sqrt(eps) 3], options);
        end
    end
    if length(beta) ~= length(O.nes)%size(O.S{O.M(1)},1)
        beta = beta(1).*ones(length(O.nes),1);%ones(size(O.S{O.M(1)},1),1);
    end
    if any(beta(O.nes) <= 0), error('All beta parameters should be positive.'); end
    O.beta = beta;                % Kotz beta parameter
    
    return
end

if strcmpi(param,'lambda')
    lambda = value;
    if isempty(lambda)
        % lambda = pi/6; % 0.523541354951965;
        % Get variance-secific lambda
        % Using the definition of alpha in the MISA paper
        lambda = zeros(length(O.nes),1);%zeros(size(O.S{O.M(1)},1),1);
        desired_stddev = pi/sqrt(3); % value that matches stddev of logistic: std = pi/sqrt(3)
        %     options = optimset('TolX',1e-100);
        for kk = find(O.nes')
            %f = @(lmb) ((lmb.^(-1./O.beta(kk)) .* gamma(O.nu(kk) + 1./O.beta(kk))) ./ (O.d(kk) .* gamma(O.nu(kk)))) - desired_stddev.^2;
            %g = @(lmb) gammaln(O.nu(kk) + 1./O.beta(kk)) - 1./O.beta(kk).*log(lmb) - log(O.d(kk)) - gammaln(O.nu(kk)) - 2*log(desired_stddev);
            %lambda(kk) = fzero(g,[0.1 100], options);
            lambda(kk) = exp(O.beta(kk) * (gammaln(O.nu(kk) + 1./O.beta(kk)) - log(O.d(kk)) - gammaln(O.nu(kk)) - 2*log(desired_stddev)));
        end
    end
    if length(lambda) ~= length(O.nes)%size(O.S{O.M(1)},1)
        lambda = lambda(1).*ones(length(O.nes),1);%ones(size(O.S{O.M(1)},1),1);
    end
    if any(lambda(O.nes) <= 0), error('All lambda parameters should be positive.'); end
    O.lambda = lambda;            % Kotz lambda parameter
    
    return
end

end