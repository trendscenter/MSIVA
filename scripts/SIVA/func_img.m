function func_sim(s_test_id)
% @gsd: generate simulated data
% @gsm: generate simulated mixing matrix

debug = false;
root_dir = "/data/users4/xli/MSIVA/MSIVA/";
addpath(genpath(root_dir));
addpath(genpath('/trdapps/linux-x86_64/matlab/toolboxes/GroupICATv4.0b/icatb/'));

res_dir = root_dir+"results/img_sz";

s_test = "s"+num2str(s_test_id)

%% load data
seed=7;
num_pc = 12;
num_iter = 0;
M = [1, 2];

% UKB data
% X = load("/data/users4/xli/MSIVA/MSIVA/results/mat/sMRI-fMRI/X_ukb.mat").X;

% patient data
X = load("/data/users4/xli/MSIVA/MSIVA/results/mat/sMRI-fMRI/X_sz.mat").X;

%% define test subspace structure
if s_test == "s1"
    St = {[1 2], [1 2]; ...
      [3 4 5], [3 4 5]; ...
      [6 7 8 9], [6 7 8 9]; ...
      [   10], [     ]; ...
      [   11], [     ]; ...
      [   12], [     ]; ...
      [     ], [   10]; ...
      [     ], [   11]; ...
      [     ], [   12]};
elseif s_test == "s2"
    St = {[1  2], [1  2]; ...
        [3  4], [3  4]; ...
        [5  6], [5  6]; ...
        [7  8], [7  8]; ...
        [9 10], [9 10]; ...
        [  11], [    ]; ...
        [  12], [    ]; ...
        [    ], [  11]; ...
        [    ], [  12]};
elseif s_test == "s3"
    St = {[1 2 3], [1 2 3]; ...
        [4 5 6], [4 5 6]; ...
        [7 8 9], [7 8 9]; ...
        [   10], [     ]; ...
        [   11], [     ]; ...
        [   12], [     ]; ...
        [     ], [   10]; ...
        [     ], [   11]; ...
        [     ], [   12]};
elseif s_test == "s4"
    St = {[1 2 3 4], [1 2 3 4]; ...
        [5 6 7 8], [5 6 7 8]; ...
        [  9], [   ]; ...
        [ 10], [   ]; ...
        [ 11], [   ]; ...
        [ 12], [   ]; ...
        [   ], [  9]; ...
        [   ], [ 10]; ...
        [   ], [ 11]; ...
        [   ], [ 12]};
elseif s_test == "s5"
    St = {[ 1], [ 1]; ...
        [ 2], [ 2]; ...
        [ 3], [ 3]; ...
        [ 4], [ 4]; ...
        [ 5], [ 5]; ...
        [ 6], [ 6]; ...
        [ 7], [ 7]; ...
        [ 8], [ 8]; ...
        [ 9], [ 9]; ...
        [10], [10]; ...
        [11], [11]; ...
        [12], [12]};
end

K = size(St,1);
eta = ones(K,1);
beta = ones(K,1);
lambda = ones(K,1);
Stest = cell(size(M));

for mm = M
    if issparse(St{mm})
        Stest{mm} = St{mm};
    else
        ii = [];
        jj = [];
        for ii_ = 1:K
            jj_ = length(St{ii_,mm});
            if jj_ ~= 0
                jj = [jj St{ii_,mm}];
                ii = [ii ii_*ones(1,jj_)];
            end
        end
        Stest{mm} = sparse(ii, jj, ones(1,sum([St{:,mm}] ~= 0)), ...
            K, sum([St{:,mm}] ~= 0), sum([St{:,mm}] ~= 0));
    end
end

%% run optimization
if s_test ~= "s5"
    % S1-S4
    num_iter = 20;
    run_optim = true;
else
    % S5/ MMIVA does not run combinatorial optimization due to a permutation error
    num_iter = 0;
    run_optim = false;
end

%% unimodal
tic
[data1_um, aux_um, isi_um, loss_um, Yinit_um, Winit_um] = run_pca_ica(X, Stest, M, num_pc, num_iter, seed, run_optim);
t_um = toc; % second

%% unimodal + multimodal
tic
[data1_ummm, aux_ummm, isi_ummm, loss_ummm, Yinit_ummm, Winit_ummm] = run_mgpca_ica(X, Stest, M, num_pc, num_iter, seed, run_optim);
t_ummm = toc; % second

%% multimodal
tic
[data1_mm, aux_mm, isi_mm, loss_mm, Yinit_mm, Winit_mm] = run_mgpca_gica(X, Stest, M, num_pc, num_iter, seed, run_optim);
t_mm = toc; % second

%% save outputs
out_dir = s_test;

W_um = data1_um.W;
W_ummm = data1_ummm.W;
W_mm = data1_mm.W;

Y_um = data1_um.Y;
Y_ummm = data1_ummm.Y;
Y_mm = data1_mm.Y;

% evaluate loss on full data
[J_um, ~] = data1_um.objective_sc_();
[J_ummm, ~] = data1_ummm.objective_sc_();
[J_mm, ~] = data1_mm.objective_sc_();

save(fullfile(res_dir,out_dir,'data1_um.mat'),'data1_um');
save(fullfile(res_dir,out_dir,'data1_ummm.mat'),'data1_ummm');
save(fullfile(res_dir,out_dir,'data1_mm.mat'),'data1_mm');
save(fullfile(res_dir,out_dir,'aux.mat'),'aux_um','aux_ummm','aux_mm');
save(fullfile(res_dir,out_dir,'loss.mat'),'loss_um','loss_ummm','loss_mm');
save(fullfile(res_dir,out_dir,'Winit.mat'),'Winit_um','Winit_ummm','Winit_mm');
save(fullfile(res_dir,out_dir,'Yinit.mat'),'Yinit_um','Yinit_ummm','Yinit_mm');
save(fullfile(res_dir,out_dir,'time.mat'),'t_um','t_ummm','t_mm');
save(fullfile(res_dir,out_dir,'Wtest.mat'),'W_um','W_ummm','W_mm');
save(fullfile(res_dir,out_dir,'Ytest.mat'),'Y_um','Y_ummm','Y_mm');
save(fullfile(res_dir,out_dir,'J.mat'),'J_um','J_ummm','J_mm');

end