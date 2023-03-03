%% Prep starting point: optimize RE to ensure initial W is in the feasible region

% Mold = data1.M;
% Sold = data1.S;
% bold = data1.beta;
% lold = data1.lambda;
% eold = data1.eta;
% old_d_k = data1.d_k;
% oldsc = data1.sc; % Save scale-control parameter
% S_ = Sold;

% Do each dataset separately
% for mm = Mold
%         % Set S to sparse identity matrix
%         S_{mm} = sparse(1:data1.C(mm),1:data1.C(mm),ones(data1.C(mm),1), ...
%         data1.C(mm),data1.C(mm),data1.C(mm)); %eye(O.C(mm));
%         % redefine parameters: accounts for changing subspace structure/size
%         ixe = cumsum(cell2mat(old_d_k(mm)));
%         ixb = [1; (ixe(1:end-1)+1)];
%         b = []; l = []; e = [];
%         for kk = 1:length(ixe)
%         b(ixb(kk):ixe(kk)) = bold(1);
%         l(ixb(kk):ixe(kk)) = lold(1);
%         e(ixb(kk):ixe(kk)) = eold(1);
%         end
%         data1.update(S_,mm,b',l',e'); % Update S and set M = mm
%         w0 = ut.stackW(W0(mm));
%         optprob = ut.getop(w0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 1e-7);
%         optprob.options.MaxIter = 1;
%         %                 data1.setREapproach('WT')
%         [woutW0,fval,exitflag,output] = fmincon(optprob);
%         optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 10^floor(log10(output.firstorderopt))/4);
%         %                 optprob.options.MaxFunEvals = 20;
%         %                 [woutW0,fval,exitflag,output] = fmincon(optprob);
%         %                 optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, optprob.options.TolX);
%         %                 %         optprob.options.MaxFunEvals = 15;
%         %                 data1.setREapproach('PINV')
%         [woutW0,fval,exitflag,output] = fmincon(optprob);
% end

% Restore original sate
% data1.update(Sold,Mold,bold,lold,eold);
woutW0 = data1.stackW(data1.W);


%% Define objective parameters and run optimization

f = @(x) data1.objective(x); % function to be optimized
% f = @(x) MISA_mean(x); 

% c = @(x) data1.con_RE(x); % nonlinear constraint function
c = [];
barr = 1; % barrier parameter
m = 1; % number of past gradients to use for LBFGS-B (m = 1 is equivalent to conjugate gradient)
N = size(X(M(1)),2); % Number of observations
Tol = .5*N*1e-9; % tolerance for stopping criteria

% Set optimization parameters and run
optprob = ut.getop(woutW0, f, c, barr, {'lbfgs' m}, Tol);
[wout,fval,exitflag,output] = fmincon(optprob);
% data1.objective(woutW0); % reset

% Prep and run combinatorial optimization
aux = {data1.W; data1.objective(ut.stackW(data1.W))};
data1.MISI(A)

% figure,imagesc(data1.W{1}*sim_siva.A{1},max(max(abs(data1.W{1}*sim_siva.A{1}))).*[-1 1]);colorbar();
% figure,imagesc(data1.W{end}*sim_siva.A{end},max(max(abs(data1.W{end}*sim_siva.A{end}))).*[-1 1]);colorbar();

data1.combinatorial_optim()
data1.MISI(A)
% for ct = 2:3
%         data1.combinatorial_optim()
%         optprob = ut.getop(ut.stackW(data1.W), f, c, barr, {'lbfgs' m}, Tol);
%         [wout,fval,exitflag,output] = fmincon(optprob);
%         aux(:,ct) = {data1.W; data1.objective_()};
%         data1.MISI(A)
% end
% [~, ix] = min([aux{2,:}]);
% data1.objective(ut.stackW(aux{1,ix}));
