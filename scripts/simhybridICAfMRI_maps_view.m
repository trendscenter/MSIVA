%% Set path to the functional data
if ispc
%     basepath = fullfile('\\loki','export','mialab');
else
    basepath = fullfile('/','home','rogers','repos');
end
datapath = fullfile(basepath,'data');
mskpath = fullfile(datapath,'masks');

% Set path to save result files
outpath = fullfile(basepath,'code','MISA','results','hyb','ICA','temporal');

addpath(genpath(fullfile('~','software','spm12')))

COMPS = 1:20;
sfile = '/home/rogers/software/MCIv4/ch2better_aligned2EPI.nii';
bname_GT = 'simhybridICAfMRI_maps_GT.nii';
fname_GT = fullfile(outpath,bname_GT);
bname_RM = 'simhybridICAfMRI_maps_RE_MISA.nii';
fname_RM = fullfile(outpath,bname_RM);

%% Load mask
Vm = spm_vol(fullfile(mskpath,'rest_hcpMask.img'));
msk = logical(spm_read_vols(Vm));
szm = size(msk);

VGT = spm_vol(fname_GT);
GT4D = spm_read_vols(VGT);
VRM = spm_vol(fname_RM);
RM4D = spm_read_vols(VRM);
sz = size(GT4D);

GT2D = reshape(GT4D,prod(sz(1:3)),sz(4));
GT2Dm = GT2D(msk(:),:)';
RM2D = reshape(RM4D,prod(sz(1:3)),sz(4));
RM2Dm = RM2D(msk(:),:)';

COMPS = 1:sz(4);

ut = utils;
RR = abs(corr(GT2Dm',RM2Dm'));
p = ut.munkres(-RR);
[~,ixs] = sort(diag(RR(:,p))','descend');

RGB = ones(256,3);
RGB(129:256,2) = linspace(1,0,128);
RGB(1:128,3) = linspace(1,0,128);
RGB(129:256,3) = zeros(128,1);
hFc = figure('Renderer', 'opengl');
imagesc(RR(ixs,p(ixs)),[0 1])
axis image; colormap(RGB); 
set(gca,'Ytick',COMPS,'Yticklabels',ixs,'Xtick',COMPS,'Xticklabels',p(ixs))
ylabel('GT')
xlabel('RM')
hcb = colorbar;
ylabel(hcb,'Correlation, r')
title('Temporal ICA of fMRI')
drawnow();

export_fig(hFc,fullfile(outpath,'simhybridICAfMRI_corr.pdf'),'-pdf','-opengl')

RM2Dms = zeros(size(RM2Dm));
for cc = COMPS
    mdl = fitlm(RM2Dm(p(cc),:)',GT2Dm(cc,:)');%,'Linear','Intercept',false);
    RM2Dms(p(cc),:) = mdl.Coefficients.Estimate('x1',:)*RM2Dm(p(cc),:);
end
out_maps_3D = to_vol(RM2Dms, msk); % 4D array

%%
fname_RMs = fullfile(outpath,'support',['k' bname_RM]);
mci_create_4DNiftifile(fname_RMs, out_maps_3D, Vm.mat)

interpnames_GT = mci_interp2struct(fname_GT, COMPS, sfile, fullfile(outpath,'support'));
interpnames_RMs = mci_interp2struct(fname_RMs, COMPS, sfile, fullfile(outpath,'support'));

%%
Vs = spm_vol(sfile);
SM = spm_read_vols(Vs);
SM = flip(SM,1); % TO make LEFT appear on the LEFT

fname_GTi = fullfile(outpath,'support',['i' bname_GT]);
Vm = spm_vol(fname_GTi);
GT = spm_read_vols(Vm);
GT = flip(GT,1); % TO make LEFT appear on the LEFT

fname_RMik = fullfile(outpath,'support',['ik' bname_RM]);
Vf = spm_vol(fname_RMik);
RM = spm_read_vols(Vf); % load in the volume
RM = flip(RM,1); % TO make LEFT appear on the LEFT
%RM(abs(RM) > 5) = 5;

load('/home/rogers/software/MCIv4/CM_coldhot_256.mat', 'CM')

sz = size(SM);
allfmaps = [];
allamaps = [];
allsmaps = [];
afrow = [];
aarow = [];
asrow = [];
cc_ = 0;
for cc = ixs
    cc_ = cc_ + 1;
    GTcc = GT(:,:,:,cc);
    RMcc = RM(:,:,:,p(cc));
    [mv,ix] = max(abs(GTcc(:)));
    sg = sign(GTcc(ix));
    [ii,jj,kk] = ind2sub(sz,ix);
    afrow = [afrow my_3views(sg*GTcc./mv,[ii,jj,kk])];
    afrow = [afrow my_3views(sg*RMcc./mv,[ii,jj,kk])];
    aarow = [aarow my_3views(abs(GTcc)>2*std(GTcc(:)),[ii,jj,kk])];
    aarow = [aarow my_3views(abs(RMcc)>2*std(RMcc(:)),[ii,jj,kk])];
    as3view = my_3views(SM,[ii,jj,kk]);
    %figure, imagesc(as3view), axis image xy; colormap gray; colorbar;
    asrow = [asrow as3view as3view];
    if mod(cc_,5) == 0
        allfmaps = [afrow; allfmaps];
        allamaps = [aarow; allamaps];
        allsmaps = [asrow; allsmaps];
        afrow = [];
        aarow = [];
        asrow = [];
    end
end

%%
hF = figure('position',[1 1 1600 821],'Renderer', 'opengl');

ax_s = axes;
imagesc(allsmaps,[0 max(SM(:))]), axis image xy
colormap(ax_s,gray)
axis(ax_s,'off')

ax_f = axes;
imagesc(allfmaps,max(abs(allfmaps(:)))*[-1 1]), axis image xy
colormap(ax_f,CM)
imh = get(ax_f,'children');
set(imh,'alphaData',double(allamaps))
axis(ax_f,'off')

%Annotations
ax_a = axes;
title('Ground-Truth (GT) vs RE+MISA Estimate')
hold(ax_a,'on')
axis(ax_a,'off')
axis image xy
total = size(allfmaps);
offy = floor(sz(2)/4);
x = .5+[0 total(2)];
y_stride = 2*sz(3) + sz(2) + 2*offy;
for ll = 1:3
    y = total(1) - .5 - y_stride*ll;
    plot(ax_a,x,[y y],'w')
end
y = .5+[0 total(1)];
x_stride = 2*sz(2);
for ll = 1:4
    x =  .5 + x_stride*ll;
    plot(ax_a,[x x],y,'w')
end
ll = 0;
offx = (sz(2)-1)/2;
x = offx;
offy2 = y_stride - offy;
r = diag(RR(:,p))';
cc_ = 0;
for ii = ixs
    cc_ = cc_ + 1;
    y = total(1) + 20.5 - offy - y_stride*(ll);
    text(ax_a,x,y,['GT' num2str(ii,'%02d')],...
        'HorizontalAlignment','center','Color','w')
    text(ax_a,x+offx,y-offy2 + 5,['r = ' num2str(r(ii),'%.4f')],...
        'HorizontalAlignment','center','Color','w')
    x = x + sz(2);
    text(ax_a,x,y,['RM' num2str(p(ii),'%02d')],...
        'HorizontalAlignment','center','Color','w')
    x = x + sz(2);
    if mod(cc_,5) == 0
        x = (sz(2)-1)/2;
        ll = ll + 1;
    end
end

linkaxes([ax_s ax_f ax_a])

cbhs = colorbar(ax_s,'eastoutside');
set(cbhs,'Visible','off')
cbhf = colorbar(ax_f,'eastoutside');
ylabel(cbhf,'arbitrary units')
cbha = colorbar(ax_a,'eastoutside');
set(cbha,'Visible','off')

export_fig(hF,fullfile(outpath,'simhybridICAfMRI_maps.pdf'),'-pdf','-opengl')

%%
% addpath('\\loki\export\mialab\users\eswar\software\utilities\')
% addpath('\\loki\export\mialab\users\eswar\software\export_fig\')
% addpath('\\loki\export\apps\linux-x86\matlab\toolboxes\Matlab_Utilities\')
% 
% [anat,~,XYZ] = icatb_read_data('\\loki\export\apps\linux-x86\matlab\toolboxes\MCIv4\ch2better_aligned2EPI.nii');
% coords = unique(XYZ(3,:));
% gsz = fliplr(coords(34:8:end-18));

% load(fullfile(outpath,[infile '.mat']),'W0','resig');


% thresh = 20;
% 
% hF = figure('color',[1 1 1],'position',[10 9 1009 988],'Name','aname');
% subplot(8,1,1);
% my_make_composite(afmap3D,anat,thresh,'',gsz,coords);
% title('atitle');
% 

%%





