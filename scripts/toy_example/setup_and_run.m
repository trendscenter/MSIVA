%% Instructions to run this script:
% % Set current folder to be the MISA folder, e.g.:
% cd('~/MISA')
% % Add the folder containing this file to the path, e.g.:
addpath('./scripts')
% % Add toy example folder to the path, e.g.:
addpath('./scripts/toy_example')
% % Clear all and RUN:
close all; clear; clc
% % run this script


% Below, orthogonal multimodal features are combined with sources.
% After setting up the subspace structure and initial unmixing matrix (W), 


%% Load the true sources (Y) and true mixing matrices (A) to generate the mixtures

% Load sources
load(fullfile('.','scripts','toy_example','jointsourcesMISA_case6.mat'), 'Sgt')
Y = Sgt;

% Load Tetris-String features
load(fullfile('.','scripts','toy_example','A_TetrisString_orthofeat_r001.mat'),'A');

% Generate mixtures
% ind = [1,3,4,5];
X{1} = A{1}*Y{1};
% X{2} = A{2}(:,ind)*Y{2}(ind,:);
X{2} = A{2}*Y{2};

%% Define subspace structure for the sources (the following matches what was used to generate the sources)

% Define the number of datasets (here, the number of modalities)
M = 1:length(X);

S = cell(1,2);           % Cell array: each cell contains a matrix K x C(m).
                         % Each k-th row has 0's and 1's to indicate what
                         % source go within the k-th subspace in dataset m

% Modality 1 = dataset 1
% S{1} = zeros(4);
% S{1}([1 6 11 15]) = ones(1,4);
S{1} = [1 0 0 0;... % source 1 into subspace 1
        0 1 0 0;... % source 2 into subspace 2
        0 0 1 0;... % source 3 into subspace 3
        0 0 0 1];   % source 4 into subspace 4

% Modality 2 = dataset 2
% S{2} = zeros(4,6);
% S{2}([1 6 10 15 20 24]) = ones(1,6);
S{2} = [1 0 0 0;... % source 1 into subspace 1
        0 1 0 0;... % source 2 into subspace 2
        0 0 1 0;... % source 3 into subspace 3
        0 0 0 1];   % source 4 into subspace 4

get_MISA_parameters


%% Initialize MISA object

data1 = MISAK(w0, M, S, X, ...
                beta, eta, lambda, ...
                gradtype, sc, preX);

% data1 = MISAKRE(w0, M, S, X, ...
%                 beta, eta, lambda, ...
%                 gradtype, sc, preX, ...
%                 REtype, REapproach, RElambda, ...
%                 REref, REreftype, REreflambda, rC);

%% Run MISA: PRE + LBFGS-B + Nonlinear Constraint + Combinatorial Optimization
execute_full_optimization

% NOTE: toggle lines 41-43 in @utils/getop.m to see detailed optimization iterations (considerably slower)


%% Check results
fprintf("\nFinal MISI: %.4f\n\n", data1.MISI(A))
% typically, a number < 0.1 indicates successful recovery of the sources


%% Visualize recovered (mixing) patterns
view_results
