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
REapproach = 'WT'; %'WT'; % 'PINV' for pseudoinverse or W

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
%     [u, s, v] = svd(randn(size(A{mm}')),'econ');
    num_comp = size(S{mm},2);
    [u, s, v] = svd(randn(size(A{mm}(:,1:num_comp)')),'econ');
    W0{mm} = u*v';
end

ut = utils;
w0 = ut.stackW(W0(M)); % vectorize unmixing matrix for compatibility with Matlab's optimization toolbox
