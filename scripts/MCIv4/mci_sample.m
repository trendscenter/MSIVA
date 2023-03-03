%% Set path to the functional data
inpath = '/export/mialab/hcp/vbm/ICA75/';

%% Set path to save interpolated files
outpath = '/export/mialab/hcp/vbm/ICA75/practice/';

%% Filename of the .nii file to make composite image of - usually directly from GIFT output
fname = 'swc1_Group_ICA75__sub01_component_ica_s1_.nii';

%% Some components of interest - here are four
COMPS = [4 50 51 75];

%% Set the full filename of structural file to use as the structural underlay
% recommended for nicest anatomical structure, and properly aligned to the EPI template that we use for normalization:
sfile = '/export/mialab/users/eallen/makeCompositeImage/ch2better_aligned2EPI.nii';

% other possible underlays:
%sfile = '/export/mialab/users/eallen/makeCompositeImage/ch2bet.nii';
%sfile = '/export/apps/linux-x86/matlab/toolboxes/spm5/canonical/single_subj_T1.nii';
%sfile = '/export/apps/linux-x86/matlab/toolboxes/spm5/templates/T1.nii';

% sfile = convertpath_linux2pc(sfile);
% inpath = convertpath_linux2pc(inpath);
% outpath = convertpath_linux2pc(outpath);

%% First, interpolate the functional data to the structural underlay dimensions -- may take a bit of time
% NOTE: in your icatb_defaults file, you must set USE_DEFAULT_SLICE_RANGE = 0;  annoying, but necessary for now

interpnames = mci_interp2struct(fullfile(inpath,fname), COMPS, sfile, outpath);


%% Now plot (in a loop) with the default settings
for ii = 1:length(interpnames)
    plot_title = ['Group Map, Componment ' num2str(COMPS(ii))];
    [slices, FH] = mci_makeimage(interpnames{ii}, sfile);
%  [slices, FH] = mci_makeimage(mapfile, structfile, ncluster, minclustersize, threshold_low, threshold_high, absflag, datalabel, units, ch_color, cmfile, slices)
%   [slices, FH] = mci_makeimage(interpnames{ii}, sfile, ncluster, minclustersize, threshold_low, [], 1, plot_title, 'intensity (a.u.)', []);
end


%% Set some parameters for the plotting
absflag = 1;
ncluster= 2;
threshold_low = 3;
minclustersize = 5000;