function [slices, FH] = mci_makeimage(mapfile, structfile, vol, varargin)
% mci_makeimage() - Overlays functional map onto structural map.  
%
% Usage:
%  >> [slices, FH] = mci_makeimage(mapfile, structfile, vol)
%  >> [slices, FH] = mci_makeimage(mapfile, structfile, vol, 'key1', val1,...)
%
% INPUTS:
% mapfile        = functional data filename (.nii file), e.g., 'imean_comps.nii'
%                  Should be in same space as structural image [see mci_inter2struct()]
% structfile     = filename of structural image (.nii file), e.g., 'chi2.nii'
% vol            = volume of functional data to overlay, e.g., 3
% 
% OPTIONAL INPUTS:
% 'ncluster'       = number of clusters at which to plot orthogonal slices [range 1-6; default = 1];
%                  Alternatively, ncluster can be set to 'i' for interactive mode where users pick the slices.
%                  In this case, the 'minclustersize' is ignored.
% 'minclustersize' = mininum size of clusters [default = 2000 (interpolated space)]
% 'threshold_low'  = value for lower threshold [default = 0.25 * max value]
% 'threshold_high' = value for high threshold [default = absolute maximum value]
% 'absflag'        = 0|1 ; 1 to threshold based on absolute magnitude and plot on symmetric colorscale [default = 1]
% 'datalabel'      = Label for each map, e.g., 'Component 1 T-map' [default = '']
% 'units'          = string for data units, e.g., 't-statistic' [default = '']
% 'ch_color'       = color for cross hairs [default is gray = [0.5, 0.5, 0.5]; 
%                   For no crosshairs, set to an empty string: ''
% 'cmfile'         = filename for colormap (.mat file), e.g.,% 'CM_green_to_red_256.mat'
%                  [default = 'CM_coldhot_256.mat'].
%                  Alternatively, cmfile can be a N x 3 matrix to be used as a colormap. N = 256 will work best.
% 'slicemethod'    = method for picking slices [default = 'dmax']
%                   'dmax': slices @ maximum point, 
%                   'dvar': slices with maximum variance,
%                   'dsum': slices with maximum sum} 
%                  Alternatively, slicemethod can be a structure with slice indices (format: S.x, S.y, S.z)
% 'slicedir'       = aspect in which to show slices; 'h' = horizontal; 'v'= vertical; [default = 'h']
% OUTPUTS:
% slices         = structure of slice indices at which images were plotted [S.x, S.y, S.z]
% FH             = figure handle
%
% See also: mci_interp2struct()

%% Load in the functional data
Vf = spm_vol(mapfile);
F = spm_read_vols(Vf(vol(1))); % load in the volume
maxval = max(abs(F(:)));

%% Load in the structural data
S = spm_read_vols(spm_vol(structfile));

%% Load in the Coordinate file
[dpath, dname] = fileparts(mapfile);
Cfname = fullfile(dpath, [dname '_coords.mat']);
load(Cfname);


%--------------------------------------------------------------------------
%% Deal with optional inputs
mci_dir = fileparts(which('mci_makeimage'));
def_cmfile = fullfile(mci_dir, 'CM_coldhot_256.mat');
%         'key'            type             range                   default
Flist = {'ncluster'  {'integer' 'string'}  {[1 6], 'i'}             1;
        'minclustersize' 'integer'         [1 numel(S)]             2000;
        'threshold_low'  'real'            [-maxval maxval]         0.25*maxval;
        'threshold_high' 'real'            [-5*maxval 5*maxval]         maxval
        'absflag'        'integer'         [0 1]                    1;
        'datalabel'      'string'          []                       '' ;
        'units'          'string'          []                       '';
        'ch_color'    {'string' 'real'}    []                       [.5 .5 .5];  
        'cmfile'      {'string' 'real'}    []                       def_cmfile;
        'slicedir'    'string'             {'h', 'v'}               'h';
        'slicemethod' {'struct' 'string'}  {'dvar', 'dmax', 'dsum'}  'dmax'};
opt = mci_finputcheck(varargin, Flist, 'mci_makeimage', 'ignore');
if ischar(opt), error(opt); end;
% all options now stored as fields in structure 'opt'
%--------------------------------------------------------------------------

%% Set any voxels greater than the max value to the max value 
F(F>opt.threshold_high) = opt.threshold_high;
fprintf('Thresholding: %0.2f to %0.2f\n', opt.threshold_low, opt.threshold_high)

%% Load the colormap
if ischar(opt.cmfile)
    load(opt.cmfile) % loads variable CM
else
    CM = opt.cmfile;
end
%% Make the images
[FH, slices] = mci_plotcomps(F, S, ... 
               opt.ncluster, opt.minclustersize, opt.threshold_low, opt.threshold_high, ... 
               opt.absflag, C, opt.datalabel, opt.units, opt.ch_color, CM, opt.slicemethod, opt.slicedir);

set(FH, 'Name', ['Volume ' num2str(vol)]);



