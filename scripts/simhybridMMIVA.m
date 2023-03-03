%function T = simtemporalIVA_mCCA(Run, part_a, part_b, inn, iaa, ir0, er0)
close all; clear all; clc
Run = 1, part_a = 1, part_b = 1, inn = 1, iaa = 1, ir0 = 1, er0 = 1
format compact
ut = utils;

if ispc
    basepath = fullfile('\\loki','export','mialab');
else
    basepath = fullfile('/','export','mialab');
end
outpath = fullfile(basepath,'users','rsilva','repos','code','MISA',...
    'results','hyb','IVA','temporal');
fname = ['temporalMMIVAhyb_w0_results_rr' num2str(Run,'%03d') ...
    '_ir0' num2str(ir0,'%03d') '_er0' num2str(er0,'%03d') '.csv'];

if exist(fullfile(outpath, fname),'file')
    T = readtable(fullfile(outpath, fname));
    ir0 = T.R0(end);
else
    T = cell2table(cell(0,12));
    T.Properties.VariableNames = {'SNRdB','Acond','Run','R0','Initial',...
        'Alg','Iter','WallTime','Cost','MISI','MMD','MMSE'};
end

K = 20; V = 20; M = 3; % K = 8, V = 8 worked with GM
% C = 1:10; V = 8000;
% S_      = [{[1 2 3];   [4]; [5 6 7];       [8]; [9]}, ...
%            {  [1 2]; [3 4];      []; [5 6 7 8]; [9]}, ...
%            {  [1 2];   [3]; [4 5 6];   [7 8 9];  []}];
% V = [2000 1000 1500];
N = 600;
Acond = [1 3 7 15];
SNR = [999   9   1   1   1; % source power
    1   1   1   9 999];% noise power
SNR = sum(SNR)./SNR(2,:); %(1+99)/99;

if part_a & part_b
    inna = inn; innb = 1;
    iaaa = iaa; iaab = 1;
    ir0a = ir0; ir0b = ir0;
    er0a = er0; er0b = er0;
else
    inna = inn; innb = inn;
    iaaa = iaa; iaab = iaa;
    ir0a = ir0; ir0b = ir0;
    er0a = er0; er0b = er0;
end
rng(100)
r0seedsa = floor(10000*rand(1,1000));
r0seedsb = floor(10000*rand(1,1000));
nn = 3;
aa = iaaa;

rr = Run;
seed = 1980+rr;

inpath{1} = fullfile(basepath,'hcp','vbm','ICA75');
in_fname{1} = 'swc1_Group_ICA75__ica.mat';
load(fullfile(inpath{1}, in_fname{1}), 'icasig')
A{1} = icasig(1:K,:)';

inpath{2} = fullfile(basepath,'hcp','groupica','rest','HC_603_S100_G75');
in_fname{2} = 'rest_hcp_ica.mat';
load(fullfile(inpath{2}, in_fname{2}), 'icasig')
A{2} = icasig(1:K,:)';

% ref_folder = '/export/mialab/users/dbridwell/fit_example/DATA/erp_fmri/ANALYSIS';
% ref_fname = 'AodfMRIEEG_ica_comb_2.mat';
% load(fullfile(ref_folder,ref_fname), 'icasig')
% A{1} = icasig(1:K,:)';

inpath{3} = fullfile(outpath,'support');
in_fname{3} = 'GICA_COBRE_30_agg__component_ica_.nii';
myV = spm_vol(fullfile(inpath{3},in_fname{3}));
D = spm_read_vols(myV);
sz = size(D);
D = reshape(D,prod(sz(1:3)),sz(4));
myVm = spm_vol(fullfile(inpath{3},'GICA_COBRE_30Mask.nii'));
msk3D = logical(spm_read_vols(myVm));
A{3} = D(msk3D(:),1:K);

sim1 = sim_hybrid_IVA(seed,K,N,A,SNR(nn),'gcp');
%sim1 = sim_basic_IVA(seed,K,V,M,N,Acond(aa),SNR(nn));
seed = rng;

[~,outfile,~] = fileparts(fname);
save(fullfile(outpath,[outfile '_GT.mat']), 'sim1')

beta = .5; gradtype = 'relative'; sc = true; RElambda = 0;
data1 = setup_basic_MISAKRE(sim1, beta, gradtype, sc, RElambda);
r0 = ir0a;

fprintf('SNR: %.3f, Acond: %.2f %.2f %.2f, Run: %d, R0: %d(%d-%d)\n', ...
    SNR(nn), cond(sim1.A{1}), cond(sim1.A{2}), cond(sim1.A{3}), ...
    rr, r0, ir0b, er0b)

W0 = cell(size(sim1.A));
for mm = sim1.M
    rng(r0seedsa(r0))
    [Q,~] = qr(randn(size(sim1.A{mm})),0);
    W0{mm} = Q';
end
data1.objective(ut.stackW(W0));

Mold = data1.M;
Sold = data1.S;
bold = data1.beta;
lold = data1.lambda;
eold = data1.eta;
old_d_k = data1.d_k;
oldsc = data1.sc; % Save scale-control parameter
S_ = Sold;
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
    optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, min(1e-8,10^floor(log10(output.firstorderopt))/4));
    %                 optprob.options.MaxFunEvals = 20;
    %                 [woutW0,fval,exitflag,output] = fmincon(optprob);
    %                 optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, optprob.options.TolX);
    %                 %         optprob.options.MaxFunEvals = 15;
    %                 data1.setREapproach('PINV')
    [woutW0,fval,exitflag,output] = fmincon(optprob);
    opt(mm).fval = fval;
    opt(mm).exitflag = exitflag;
    opt(mm).output = output;
    data1.RE_()
    output
end
data1.update(Sold,Mold,bold,lold,eold);
woutW0 = data1.stackW(data1.W);
data1.objective(woutW0);

% Things to save:
W0re = data1.W;
resig = data1.Y;
step1 = data1;

% Erase
tmpX = step1.X;
tmpY = step1.Y;
step1.X = {};
step1.Y = {};

[~,outfile,~] = fileparts(fname);
save(fullfile(outpath,[outfile '_RE.mat']), ...
        'W0re','resig','step1','opt','inpath','in_fname','outpath')
step1.X = tmpX;
step1.Y = tmpY;

tic
data2 = MISAKRE(ut.stackW(cellfun(@(c) eye(c), num2cell(data1.C), 'Un', 0)), ...
    data1.M, data1.S, data1.Y, data1.beta, data1.eta, data1.lambda, data1.gradtype, data1.sc, data1.preX, ...
    data1.REtype, data1.REapproach, data1.RElambda, ...
    data1.REref, data1.REreftype, data1.REreflambda, data1.rC);
barrier = 1; m = 10;
optprob = ut.getop(ut.stackW(data2.W), @(x) data2.objective(x), [], ...
    barrier, {'lbfgs' m}, 1e-9);
% optprob.options.MaxIter = 1000;
[wout,fval,exitflag,output] = fmincon(optprob);
clear opt
opt.fval = fval;
opt.exitflag = exitflag;
opt.output = output;
output

% [w, fval, iters] = myline_search_batch(ut.stackW(data2.W),data2);
t = toc;
%         figure, imagesc(corr(data2.Y{3}',sim1.Y{3}'), [-1 1])
%         figure, imagesc(corr(data2.Y{2}',sim1.Y{2}'), [-1 1])
%         figure, imagesc(corr(data2.Y{1}',sim1.Y{1}'), [-1 1])
%         ut.MISI(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0), sim1.A, data1.S)
data2.updatesc(false)
data2.objective(ut.stackW(data2.W));
data2.updatesc(true)
data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
T2 = cell2table({10*log10(SNR(nn)), max(sim1.Acond), rr, r0, {'W0'}, ...
    {'RE+MISA'}, output.iterations, t, data2.objective_(), ...
    data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
T2.Properties.VariableNames = T.Properties.VariableNames;
T = [T;T2];
fprintf('MISA: %.5f\n', data1.MISI(sim1.A))

% Things to save:
step2 = data2;
Wica = step2.W;
Aica = cellfun(@(x) inv(x), step2.W, 'Un', 0);
W = cellfun(@mtimes, step2.W, W0re, 'Un', 0);
A = cellfun(@(x) ((x*x')\x)', W, 'Un', 0);
Ah = cellfun(@(y,x) ((y*y')\(y*x'))', step2.Y, step1.X, 'Un', 0);
icasig = step2.Y;

% Erase data from MISA object
tmpX = step2.X;
tmpY = step2.Y;
step2.X = {};
step2.Y = {};

[~,outfile,~] = fileparts(fname);
save(fullfile(outpath,[outfile '_MISA.mat']), ...
        'Wica','Aica','W0re','W','A','Ah','icasig','step2','opt','outpath')
step2.X = tmpX;
step2.Y = tmpY;

writetable(T, fullfile(outpath, fname))

% Save maps, modality 1: ------------------------------------------------
% Load mask
Vm = spm_vol(fullfile(inpath{1},'swc1_Group_ICA75_Mask.img'));
msk = logical(spm_read_vols(Vm));

% GT maps
out_maps_3D = to_vol(sim1.A{1}', msk); % 4D array
fname = fullfile(outpath,'simhybridMMIVA_maps_m1_GT.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)

% Estimated maps
out_maps_3D = to_vol(A{1}', msk); % 4D array
fname = fullfile(outpath,'simhybridMMIVA_maps_m1_RE_MISA.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)

% Save maps, modality 2: ------------------------------------------------
% Load mask
Vm = spm_vol(fullfile(inpath{2},'rest_hcpMask.img'));
msk = logical(spm_read_vols(Vm));

% GT maps
out_maps_3D = to_vol(sim1.A{2}', msk); % 4D array
fname = fullfile(outpath,'simhybridMMIVA_maps_m2_GT.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)

% Estimated maps
out_maps_3D = to_vol(A{2}', msk); % 4D array
fname = fullfile(outpath,'simhybridMMIVA_maps_m2_RE_MISA.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)

% Save maps, modality 3: ------------------------------------------------
% Load mask
Vm = spm_vol(fullfile(inpath{3},'GICA_COBRE_30Mask.nii'));
msk = logical(spm_read_vols(Vm));

% GT maps
out_maps_3D = to_vol(sim1.A{3}', msk); % 4D array
fname = fullfile(outpath,'simhybridMMIVA_maps_m3_GT.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)

% Estimated maps
out_maps_3D = to_vol(A{3}', msk); % 4D array
fname = fullfile(outpath,'simhybridMMIVA_maps_m3_RE_MISA.nii');
mci_create_4DNiftifile(fname, out_maps_3D, Vm.mat)

% data2 = MISAKRE(ut.stackW(cellfun(@(c) eye(c), num2cell(data1.C), 'Un', 0)), ...
%     data1.M, data1.S, data1.Y, data1.beta, data1.eta, data1.lambda, data1.gradtype, data1.sc, data1.preX, ...
%     data1.REtype, data1.REapproach, data1.RElambda, ...
%     data1.REref, data1.REreftype, data1.REreflambda, data1.rC);
% data2.update(Sold,Mold,bold*2,lold./2,eold);% Lap --> Gauss
% 
% %data1.objective(woutW0);
% W_ = mcca_genvar(cat(3,data2.X{:}),data2.C(1));
% for mm = data2.M
%     W{mm} = W_(:,:,mm);
% end
% data2.objective(ut.stackW(W))
% data2.updatesc(false)
% data2.objective(ut.stackW(data2.W));
% data2.updatesc(true)
% data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
% T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
%     {'RE+MISA'}, output.iterations, t, data2.objective_(), ...
%     data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
% T2.Properties.VariableNames = T.Properties.VariableNames;
% T = [T;T2];
% fprintf('mCCA: %.5f\n', data1.MISI(sim1.A))
% 
% Mold = data2.M;
% Sold = data2.S;
% bold = data2.beta;
% lold = data2.lambda;
% eold = data2.eta;
% old_d_k = data2.d_k;
% oldsc = data2.sc; % Save scale-control parameter
% S_ = Sold;
% for mm = Mold
%     % Set S to sparse identity matrix
%     S_{mm} = sparse(1:data2.C(mm),1:data2.C(mm),ones(data2.C(mm),1), ...
%         data2.C(mm),data2.C(mm),data2.C(mm)); %eye(O.C(mm));
%     % redefine parameters: accounts for changing subspace structure/size
%     ixe = cumsum(cell2mat(old_d_k(mm)));
%     ixb = [1; (ixe(1:end-1)+1)];
%     b = []; l = []; e = [];
%     for kk = 1:length(ixe)
%         b(ixb(kk):ixe(kk)) = bold(1);
%         l(ixb(kk):ixe(kk)) = lold(1);
%         e(ixb(kk):ixe(kk)) = eold(1);
%     end
%     data2.update(S_,mm,b',l',e'); % Update S and set M = mm
%     
%     barrier = 1; m = 10;
%     optprob = ut.getop(ut.stackW(data2.W(mm)), @(x) data2.objective(x), [], barrier, {'lbfgs' m}, 1e-10);
%     [wout,fval,exitflag,output] = fmincon(optprob);
% 
% end
% data2.update(Sold,Mold,bold,lold,eold);
% woutW0 = data2.stackW(data2.W);
% data2.updatesc(false)
% data2.objective(ut.stackW(data2.W));
% data2.updatesc(true)
% data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
% T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
%     {'RE+MISA'}, output.iterations, t, data2.objective_(), ...
%     data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
% T2.Properties.VariableNames = T.Properties.VariableNames;
% T = [T;T2];
% fprintf('MISA: %.5f\n', data1.MISI(sim1.A))