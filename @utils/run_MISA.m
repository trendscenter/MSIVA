function [wout,fval,exitflag,output] = run_MISA(misa_obj,W0)

%% Define objective parameters and run optimization
f = @(x) misa_obj.objective(x); % function to be optimized
% c = @(x) misa_obj.con_RE(x); % nonlinear constraint function
barr = 1; % barrier parameter
m = 1; % number of past gradients to use for LBFGS-B (m = 1 is equivalent to conjugate gradient)
N = misa_obj.N; % Number of observations
Tol = .5*N*1e-9; % tolerance for stopping criteria

ut = utils;
w0 = ut.stackW(W0); % vectorize unmixing matrix for compatibility with Matlab's optimization toolbox

% Set optimization parameters and run
optprob = ut.getop(w0, f, [], barr, {'lbfgs' m}, Tol); % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< wout !!!
% optprob.options.Display = 'iter-detailed';
[wout,fval,exitflag,output] = fmincon(optprob);
misa_obj.objective(wout); % FINAL