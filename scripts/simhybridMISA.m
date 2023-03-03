%% Instructions to run this script:
% % Set current folder to be the MISA folder, e.g.:
% cd('~/MISA')
% % Add the folder containing this file to the path, e.g.:
% addpath('./scripts')
% % Clear all and RUN:
% close all; clear all; clc
% call this script


% Below, orthogonal multimodal features are combined with sources.
% After setting up the subspace structure and initial unmixing matrix (W), 


%% Load the true sources (Y) and true mixing matrices (A) to generate the mixtures

% Load sources
load(fullfile('..','MISA-data','Sgt','case6_','jointsourcesMISA_case6.mat'), 'Sgt')
Y = Sgt;

% ERP features
load(fullfile('..','MISA-data','real_data_features','orthofeat_ERP.mat'), 'ERP')
A{1} = ERP;

% fMRI features
load(fullfile('..','MISA-data','real_data_features','orthofeat_fMRI.mat'), 'fMRI')
A{2} = fMRI;

% Generate mixtures
X{1} = A{1}*Y{1};
X{2} = A{2}*Y{2};


%% Define subspace structure for the sources (the following matches what was used to generate the sources)

% Define the number of datasets (here, the number of modalities)
M = 1:length(X);

S = cell(1,2);           % Cell array: each cell contains a matrix K x C(m).
                         % Each k-th row has 0's and 1's to indicate what
                         % source go within the k-th subspace in dataset m

% Modality 1 = dataset 1
%S{1} = zeros(4);
%S{1}([1 6 11 15]) = ones(1,4);
S{1} = [1 0 0 0;... % source 1 into subspace 1
        0 1 0 0;... % source 2 into subspace 2
        0 0 1 1;... % sources 3 and 4 into subspace 3
        0 0 0 0];   % no sources from modality 1 into subspace 4

% Modality 2 = dataset 2
%S{2} = zeros(4,6);
%S{2}([1 6 10 15 20 24]) = ones(1,6);
S{2} = [1 0 0 0 0 0;... % source 1 into subspace 1
        0 1 1 0 0 0;... % sources 2 and 3 into subspace 2
        0 0 0 1 0 0;... % source 4 into subspace 3
        0 0 0 0 1 1];   % sources 5 and 6 into into subspace 4 

K = size(S{M(1)},1);   % Number of subspaces

% Set Kotz parameters to multivariate lapace
eta = ones(K,1);
beta = .5*ones(K,1);
lambda = ones(K,1);


%% Set additional parameters for MISA

% Use relative gradient
gradtype = 'relative';

% Enable scale control
sc = 1;

% Turn off preprocessing (still removes the mean of the data)
preX = false;


%% Set reconstruction error options:

% Use normalized MSE
REtype = 'NMSE';

% Use the transpose of W as the reconstruction approach
REapproach = 'WT'; % 'PINV' for pseudoinverse or W

% Tolerance level (0 means least error possible)
RElambda = 5e-3;

% Other parameters required by the @MISAKRE API but not used
REref = {};
REreftype = 'linearreg';
REreflambda = {.9};
rC = {[],[]};


%% Define the starting point (for this problem an orthogonal unmixing matrix, since the features A are orthogonal).

rng(100) % set the seed for reproducibility
W0 = cell(size(A));
for mm = M
    [u, s, v] = svd(randn(size(A{mm}')),'econ');
    W0{mm} = u*v';
end

ut = utils;
w0 = ut.stackW(W0(M)); % vectorize unmixing matrix for compatibility with Matlab's optimization toolbox


%% Initialize MISA object

data1 = MISAKRE(w0, M, S, X, ...
                beta, eta, lambda, ...
                gradtype, sc, preX, ...
                REtype, REapproach, RElambda, ...
                REref, REreftype, REreflambda, rC);


%% Prep starting point: optimize RE to ensure initial W is in the feasible region

Mold = data1.M;
Sold = data1.S;
bold = data1.beta;
lold = data1.lambda;
eold = data1.eta;
old_d_k = data1.d_k;
oldsc = data1.sc; % Save scale-control parameter
S_ = Sold;

% Do each dataset separately
for mm = Mold
        % Set S to sparse identity matrix
        S_{mm} = sparse(1:data1.C(mm),1:data1.C(mm),ones(data1.C(mm),1), ...
        data1.C(mm),data1.C(mm),data1.C(mm)); %eye(O.C(mm));
        % redefine parameters: accounts for changing subspace structure/size
        ixe = cumsum(cell2mat(old_d_k(mm)));
        ixb = [1; (ixe(1:end-1)+1)];
        b = []; l = []; e = [];
        for kk = 1:length(ixe)
        b(ixb(kk):ixe(kk)) = bold(1);
        l(ixb(kk):ixe(kk)) = lold(1);
        e(ixb(kk):ixe(kk)) = eold(1);
        end
        data1.update(S_,mm,b',l',e'); % Update S and set M = mm
        w0 = ut.stackW(W0(mm));
        optprob = ut.getop(w0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 1e-7);
        optprob.options.MaxIter = 1;
        %                 data1.setREapproach('WT')
        [woutW0,fval,exitflag,output] = fmincon(optprob);
        optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 10^floor(log10(output.firstorderopt))/4);
        %                 optprob.options.MaxFunEvals = 20;
        %                 [woutW0,fval,exitflag,output] = fmincon(optprob);
        %                 optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, optprob.options.TolX);
        %                 %         optprob.options.MaxFunEvals = 15;
        %                 data1.setREapproach('PINV')
        [woutW0,fval,exitflag,output] = fmincon(optprob);
end
% Restore original sate
data1.update(Sold,Mold,bold,lold,eold);
woutW0 = data1.stackW(data1.W);


%% Define objective parameters and run optimization

f = @(x) data1.objective(x); % function to be optimized
c = @(x) data1.con_RE(x); % nonlinear constraint function
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
for ct = 2:3
        data1.combinatorial_optim()
        optprob = ut.getop(ut.stackW(data1.W), f, c, barr, {'lbfgs' m}, Tol);
        [wout,fval,exitflag,output] = fmincon(optprob);
        aux(:,ct) = {data1.W; data1.objective_()};
        data1.MISI(A)
end
[~, ix] = min([aux{2,:}]);
data1.objective(ut.stackW(aux{1,ix}));

%% Correct sources within subspaces:
%% Run ICA on subspace-specific columns of estimated A{mm}
for mm = data1.M
    for kk = find(data1.nes)'
        d_k = data1.d_k{mm}(kk);
        if d_k > 1
            w0_ = ut.stackW({eye(d_k)});
            M_ = 1;
            Stmp = {sparse(1:d_k,1:d_k,ones(1,d_k),d_k,d_k)};
            X_ = {data1.W{mm}(find(full(data1.S{mm}(kk,:))),:)};
            detX_ = prod(eig((cov(X_{M(1)}'))));
            tmp = MISAK(w0_, M_, Stmp, X_, ...
                data1.beta(kk)*ones(d_k,1), data1.eta(kk)*ones(d_k,1), data1.lambda(kk)*ones(d_k,1), ...
                data1.gradtype, data1.sc, []);
            
            f = @(x) tmp.objective(x); % function to be optimized
            c = []; % nonlinear constraint function
            barr = 1; % barrier parameter
            m = 1; % number of past gradients to use for LBFGS-B (m = 1 is equivalent to conjugate gradient)
            N = size(X_(M(1)),2); % Number of observations
            Tol = 1e-11; % tolerance for stopping criteria

            % Set optimization parameters and run
            optprob_ = ut.getop(w0_, f, c, barr, {'lbfgs' m}, Tol);
            [wout_,fval_,exitflag_,output_] = fmincon(optprob_);

            detY_ = prod(eig((cov(tmp.Y{M(1)}'))));
            tmp.objective(((sqrt(detX_)/sqrt(detY_))^(1/d_k))*wout_);
            data1.W{mm}(find(full(data1.S{mm}(kk,:))),:) = tmp.Y{M(1)};
        else
            kk
            disp('No')
        end
    end
end
data1.objective(ut.stackW(data1.W))

% Save results
T = cell2table({data1.objective_(), ...
    data1.MISI(A), data1.MMD(A), data1.MMSE(Y) });
T.Properties.VariableNames = {'Cost','MISI','MMD','MMSE'};

% Things to save:
step1 = data1;
W = step1.W;
A_tr = cellfun(@(w) w', W, 'Un', 0);
A_pinv = cellfun(@(w) ((w*w')\w)', W, 'Un', 0);
A_ls = cellfun(@(y,x) ((y*y')\(y*x'))', step1.Y, step1.X, 'Un', 0);
icasig = step1.Y;

% Erase data from MISA object
tmpX = step1.X;
tmpY = step1.Y;
step1.X = {};
step1.Y = {};

%% Set input paths
basepath = '.';
scripath = fullfile(basepath,'scripts');
mskpath = fullfile('..','MISA-data','real_data_features');
% Set output path
outpath = fullfile(basepath,'results','hyb_','MISA');

fname = 'MMMISAhyb_w0_results';
save(fullfile(outpath,[fname '_MISA.mat']), ...
        'W','A','A_tr','A_pinv','A_ls','icasig','step1','outpath')

% Restore data to MISA object
step1.X = tmpX;
step1.Y = tmpY;

writetable(T, fullfile(outpath, [fname '_MISA.csv']))

addpath ~/GSU/software/spm12/
addpath('./scripts/MCIv4')
% Save maps, modality 2: ------------------------------------------------
% Load mask
Vm = spm_vol(fullfile(mskpath,'rest_hcpMask.img'));
msk = logical(spm_read_vols(Vm));

% GT maps
out_maps_3D = to_vol(A{2}', msk); % 4D array
fname = fullfile(outpath,'simhybridMISA_maps_m2_GT.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)

% Estimated maps
out_maps_3D = to_vol(A_tr{2}', msk); % 4D array
fname = fullfile(outpath,'simhybridMISA_maps_m2_RE_WT_MISA.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)


% View ERPs side by side:
% GT ERP:
figure
mm = 1;
for ii = 1:data1.C(mm)
    subplot(data1.C(mm),2,2*ii-1)
    plot(A{mm}(:,ii)), axis tight
    if ii == 1, title('Ground-Truth'); end
    ylabel(['ERP ' num2str(ii)])
    subplot(data1.C(mm),2,2*ii)
    plot(data1.W{mm}(ii,:)), axis tight
    if ii == 1, title('Estimate'); end
end

export_fig(gcf,fullfile(outpath,'simhybridMISA_timecourses.pdf'),'-pdf','-opengl')

simhybridMISA_maps_view

