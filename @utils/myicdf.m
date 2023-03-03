function X = myicdf(name, U, varargin)

X = zeros(size(U));

switch(lower(name))
    case 'double_logistic'
        neU = numel(U);
        X = zeros(size(U));
%         X = varargin{3}*my_icdf('logistic',U,varargin{1:2})+(1-varargin{3})*my_icdf('logistic',U,varargin{4:5});
        for kk = 1:neU
            X(kk) = find_icdf_dlog(U(kk),varargin{:});
        end
    case 'double_cos'
        neU = numel(U);
        X = zeros(size(U));
        for kk = 1:neU
            X(kk) = find_icdf_dcos(U(kk));
        end
    case 'normal'
        m = varargin{1};
        s = varargin{2};
        X = icdf('norm',U,m,s);
    case 'uniform'
        a = varargin{1};
        b = varargin{2};
        X = (b-a)*(U - .5) + mean([a,b]);
    case 'gev'
        e = varargin{1};
        s = varargin{2};
        m = varargin{3};
        X = icdf('gev',U,e,s,m);
    case 'scale_t'
        t = varargin{1};
        s = varargin{2};
        X = s*icdf('t',U,t);
    case 'bimodal_scale_t'
        neU = numel(U); %k = round(neU*(14/15));
        X = zeros(size(U));
        for kk = 1:neU
            X(kk) = find_icdft(U(kk),varargin{:});
        end
    case 'bimodal_gamma'
        neU = numel(U);
        X = zeros(size(U));
        for kk = 1:neU
            X(kk) = find_icdfgam(U(kk),varargin{:});
        end
    case 'gn'
        neU = numel(U);
        X = zeros(size(U));
        for kk = 1:neU
            X(kk) = find_icdfgn(U(kk),varargin{:});
        end
    case 'laplace'
        U = U - .5;
        m = varargin{1};
        sd = varargin{2};
        if sd <= 0
            error('Parameter value out of valid range.');
        end
        X = m - sd*sign(U).*log(1 - 2*abs(U));
    case 'loglap'
        U = ((.98.*(U-1/2))+1/2);
        a = varargin{1};
        b = varargin{2};
        d = varargin{3};
        q = a/(a+b);
        X(U <= q) = d*((U(U <= q).*((a+b)/a)).^(1/b));
        X(U > q) = d*((((U(U > q)-1).*(-1)).*((a+b)/b)).^(-1/a));
    case 'logistic'
        m = varargin{1};
        b = varargin{2};
        if b <= 0
            error('Parameter value out of valid range.');
        end
        X = m + b*log( U./(1-U) );
    case 'hypbsec'
        X = (-2/pi)*asinh( cot(pi*U) );
    case 'cauchy'
        m = varargin{1};
        half_iqr = varargin{2};
        X = m +half_iqr*tan( pi*(U-.5) ); %Rescaled uniform variates to [-.98*.5, .98*.5]
    otherwise
        clear X
        error('Invalid distribution name.');
end

end

function out = find_icdft(u,t1,s1,t2,s2,m2)
fun = @(x) ((14/15)*cdf('t',x/s1,t1) + (1/15)*cdf('t',(x-m2)./s2,t2)) - u;
x0 = s1*icdf('t',u,t1);
options = optimset('FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',10000, 'TolX',1e-25);
out = fzero(fun,x0,options);
end

function out = find_icdfgam(u,a,b,m)
fun = @(x) .5*cdf('gam',x+(a-1)*b-m,a,b) + .5*(1-cdf('gam',-x+(a-1)*b-m,a,b)) - u;
x0 = 0;
%if icdf('gam',u,a,b)-(a-1)*b+m;
options = optimset('FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',10000, 'TolX',1e-25);
out = fzero(fun,x0,options);
end

function out = find_icdfgn(u,m,a,b)
fun = @(x) gncdf(x,m,a,b) - u;
x0 = m;
options = optimset('FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',10000, 'TolX',1e-25);
out = fzero(fun,x0,options);
end

function out = find_icdf_dcos(u)
fun = @(x) sin(4*pi*x + pi)/(4*pi) + x - u;
x0 = .5;
options = optimset('FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',10000, 'TolX',1e-25);
out = fzero(fun,x0,options);
end

function out = find_icdf_dlog(u,m1,s1,p1,m2,s2,p2)
fun = @(x) p1.*(1./(1+exp(-(x-m1)./s1))) + p2.*(1./(1+exp(-(x-m2)./s2))) - u;
x0 = 0;
%if icdf('gam',u,a,b)-(a-1)*b+m;
options = optimset('FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',10000, 'TolX',1e-25);
out = fzero(fun,x0,options);
end