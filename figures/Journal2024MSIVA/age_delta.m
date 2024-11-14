% Load functions
addpath("/data/users4/xli/MSIVA/FSLNets");

% Load variables
ver = 4;
A = load("/data/users4/xli/MSIVA/mat/A_ukb.mat").A;
A_msiva_s2 = squeeze(A(2,2,:,:,:));
Y_msiva_s2 = load("/data/users4/xli/MSIVA/mat/ummm_neuroimaging_Y.mat").Y2;

[n_modality, ~, n_subject] = size(Y_msiva_s2);
[~, n_voxel, ~] = size(A_msiva_s2);

ay = {5*n_modality + 4}; % 5 shared subspaces x 2 modalities + 2 unique subspaces x 2 modalities

for s = 1:5
    for m = 1:2
        a = squeeze(A_msiva_s2(m,:,2*s-1:2*s));
        y = squeeze(Y_msiva_s2{m}(2*s-1:2*s,:));
        ay{2*(s-1)+m} = a * y;
    end
end

for m = 1:2
    a = squeeze(A_msiva_s2(m,:,:));
    y = squeeze(Y_msiva_s2{m});
    ay{10+2*m-1} = a(:,11) * y(11,:);
    ay{10+2*m} = a(:,12) * y(12,:);
end

AY = cat(3,ay{:});

load("/data/users4/xli/MSIVA/mat/UKB_MMIVA_C30_preregSite_SMRI_MancovanOuts_wX_FINAL.mat");
load("/data/users4/xli/MSIVA/mat/UKB_MMIVA_C30_preregSite_ALFF_MancovanOuts_wX_FINAL.mat");
load("/data/users4/xli/MSIVA/mat/UKB_MMIVA_C30_preregSite_DMRI_MancovanOuts_wX_FINAL.mat");

ageCol=find(strcmp('age_when_attended_assessment_centre_f21003_2_0',MODELUKB0s_ful.names));
age = MODELUKB0s_ful.X(:,ageCol);

sexCol = find(strcmp('sex_f31_0_0',MODELUKB0s_ful.names));
sex = MODELUKB0s_ful.X(:,sexCol);

ssnCol=find(strcmp('rSpatNorm' ,MULT0s_ful.final_terms));
fsnCol=find(strcmp('rSpatNorm' ,MULT1s_ful.final_terms));
dsnCol = find(strcmp('rSpatNorm' ,MULT2s_ful.final_terms));
ssn = MODELUKB0s_ful.X(:,ssnCol);
fsn = MODELUKB1s_ful.X(:,fsnCol);
dsn = MODELUKB2s_ful.X(:,dsnCol);

fdCol=find(strcmp('meanFD', MULT1s_ful.final_terms));
fd = MODELUKB1s_ful.X(:,fdCol);

%%
Y = age;
y = nets_demean(Y);
yy = nets_demean([y  nets_unconfound(nets_demean(y.^2),y)]); 
yy = (yy * std(yy(:,1))) ./ std(yy);
yyy = nets_demean([yy nets_unconfound(nets_demean(y.^3),yy)]); 
yyy = (yyy * std(yyy(:,1))) ./ std(yyy);
sex_age = sex .* yyy;
confY = cat(2, yyy, sex, sex_age, fd, ssn, fsn);

fpvalue = zeros(n_voxel, 1);
tpvalue = zeros(n_voxel, 14);
delta2_ = {n_voxel}; %zeros(n_voxel, n_subject, 14);
delta2p_ = {n_voxel}; %zeros(n_voxel, n_subject, 14);
beta1_ = {n_voxel};
predictor_ = {n_voxel};

%%
for ind_vox = 1:n_voxel
    
    disp(ind_vox)
    
    X = squeeze(AY(ind_vox,:,:));
    X_svd_norm = zeros([length(X), 5]);
    
    for s = 1:5 % per subspace
        X_subspace = X(:,[2*s-1,2*s]); % subject x two modalities
        X_norm = nets_normalise(X_subspace,1);
        C = X_norm' * X_norm;
        [U, S] = eig(C);
        X_svd = X_norm * U(:,2); % sign ambiguity, make it positively associated with age
        X_svd_sign = X_svd * sign( corr(X_svd, age) );
        X_svd_norm(:,s) = nets_normalise(X_svd_sign, 1);
    end
    
    % Compute residuals
    X_res_m1 = nets_nodepartial(cat(2, X_svd_norm, X(:,[1:2:9,11:12]))); % replace X with SVD+X, extract SCVs to remove SVD from X unique subspace
    X_res_norm_m1 = nets_normalise(X_res_m1);
    X_res_m2 = nets_nodepartial(cat(2, X_svd_norm, X(:,[2:2:10,13:14]))); % replace X with SVD+X, extract SCVs to remove SVD from X unique subspace
    X_res_norm_m2 = nets_normalise(X_res_m2);
    
    pcaicaU=cat(2, X_svd_norm, X_res_norm_m1(:,6:end), X_res_norm_m2(:,[11,12])); % sMRI
    
    % Normalize and partialize expression levels
    pcaicaU=nets_normalise(pcaicaU);
    pcaicaUP=nets_normalise(nets_nodepartial(pcaicaU));
    
    X=pcaicaU;
    partialX=pcaicaUP;
    predictor_{ind_vox}=partialX;

    beta1=pinv(partialX)*y;
    model=fitlm(partialX,y);
    beta1=model.Coefficients.Estimate(2:end);
    beta1_{ind_vox}=beta1;
    tpvalue(ind_vox,:)=model.Coefficients.pValue(2:end)';
    fpvalue(ind_vox)=model.ModelFitVsNullModel.Pvalue;
    
    Yb1=partialX.*beta1'; % delta 1
    delta2ica=nets_unconfound(Yb1,confY); % no need to subtract y from Yb1, given last step
    delta2icaP=nets_nodepartial(delta2ica); % partialise deltas
    delta2_{ind_vox}=delta2ica;
    delta2p_{ind_vox}=delta2icaP;
end

predictor = cat(3,predictor_{:});
delta2 = cat(3,delta2_{:});
delta2p = cat(3,delta2p_{:});
beta1 = cat(2,beta1_{:});

save(['mat/v',num2str(ver),'/delta2p_raw/predictor_sMRI.mat'], 'predictor', '-v7.3');
save(['mat/v',num2str(ver),'/delta2p/tpvalue_sMRI.mat'], 'tpvalue', '-v7.3');
save(['mat/v',num2str(ver),'/delta2p_raw/delta2_sMRI.mat'], 'delta2', '-v7.3');
save(['mat/v',num2str(ver),'/delta2p_raw/delta2p_sMRI.mat'], 'delta2p', '-v7.3');
save(['mat/v',num2str(ver),'/delta2p/beta1_sMRI.mat'], 'beta1', '-v7.3');

for s = 1:14 % per predictor
    delta2p_subspace = squeeze(delta2p(:,s,:));
    delta2p_mean = squeeze(mean(delta2p_subspace, 1));
    delta2p_median = squeeze(median(delta2p_subspace, 1));
    delta2p_std = squeeze(std(delta2p_subspace, [], 1));
    save(['mat/v',num2str(ver),'/delta2p_raw/delta2p_predictor',num2str(s),'_sMRI.mat'], 'delta2p', '-v7.3');
    save(['mat/v',num2str(ver),'/delta2p/delta2p_mean_predictor',num2str(s),'_sMRI.mat'], 'delta2p_mean', '-v7.3');
    save(['mat/v',num2str(ver),'/delta2p/delta2p_median_predictor',num2str(s),'_sMRI.mat'], 'delta2p_median', '-v7.3');
    save(['mat/v',num2str(ver),'/delta2p/delta2p_std_predictor',num2str(s),'_sMRI.mat'], 'delta2p_std', '-v7.3');
end