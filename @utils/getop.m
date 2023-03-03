function optprob = getop(w0,varargin)

load('optprob4.mat')
if ~isempty(varargin)
    if ~isempty(varargin{1})
        optprob.objective = varargin{1};
    end
    if length(varargin) > 1 && ~isempty(varargin{2})
        optprob.nonlcon = varargin{2};
    end
    if length(varargin) > 2 && ~isempty(varargin{3})
        optprob.options.InitBarrierParam = varargin{3};
    else
        optprob.options.InitBarrierParam = 0.1;
    end
    if length(varargin) > 3 && ~isempty(varargin{4})
        optprob.options.Hessian = varargin{4};
    else
        optprob.options.Hessian = {'lbfgs' 15};
    end
    if length(varargin) > 4 && ~isempty(varargin{5})
        optprob.options.TolFun = varargin{5};
        optprob.options.TolX = varargin{5}; %sqrt(eps);%5e-17;
    else
        optprob.options.TolFun = 1e-4;%1.1e-3; %sqrt(eps);%5.5e-4;
        optprob.options.TolX = 1e-9;%sqrt(eps); %sqrt(eps);%5e-17;
    end
    if length(varargin) > 5 && ~isempty(varargin{6})
        myplot = varargin{6};
    else
        myplot = @optimplotfval;
    end
end
optprob.x0 = w0(:);
optprob.lb = -100*ones(size(w0));
optprob.ub = 100*ones(size(w0));
optprob.options.TypicalX = .1*ones(size(w0));
optprob.options.MaxFunEvals = 50000;
optprob.options.MaxIter = 150;
% optprob.options.ObjectiveLimit = 0;
% optprob.options.Diagnostics = 'on';
% optprob.options.Display = 'iter-detailed';
% optprob.options.PlotFcns = {@optimplotconstrviolation @optimplotstepsize @optimplotfirstorderopt};
% optprob.options.PlotFcns = {myplot @optimplotconstrviolation @optimplotstepsize @optimplotfirstorderopt};


end