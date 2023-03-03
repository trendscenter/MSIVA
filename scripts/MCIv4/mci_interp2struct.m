function [fname, Cfname] = mci_interp2struct(MAPfile, volumes, structFile, out_path)
% mci_interp2struct() - Interpolates the functional data to the structural (underlay) dimensions.
%  
% Usage: >> [fname, Cfname] = mci_interp2struct(MAPfile, volumes, structFile, out_path);
%   
% INPUTS:
% MAPfile    = Functional data filename, should be .nii 
% volumes    = volumes of MAPfile to interpolate (typically the components of interest)
% structFile = Structural data filename, also .nii 
% out_path   = path to save interpolated files
%
% OUTPUTS:
% fname      = filenames of the interpolated maps.  Same name as input file with an 'i' prefix.
% Cfname     = filename of the interpolated map coordinates.
%
% See also: mci_makeimage()

[dpath, dname, dext] = fileparts(MAPfile);

if nargin < 4
    out_path = dpath;
end
       
%% load in the structural info
[structuralImage, structHInfo, XYZ] = icatb_read_data(structFile);
clear structuralImage
XYZ = reshape(XYZ,[3,structHInfo.DIM]);


%% interpolate
Fimages = mci_interpdata(structHInfo.V(1), MAPfile, volumes);
Fimages = permute(Fimages, [2 3 4 1]);

%% Save the coordinates in .mat format
C = mci_simplify_coords(XYZ);
Cfname = fullfile(out_path, ['i' dname '_coords.mat']);
fprintf('Saving interpolated coordinates to %s...\n', Cfname)
save(Cfname, 'C');

%% Save the interpolated data in nifti format
fname = fullfile(out_path, ['i' dname '.nii']);
fprintf('Saving interpolated data to %s...\n\t new dimensions [%d x %d x %d]...', fname, size(Fimages,1), size(Fimages,2), size(Fimages,3))
mci_create_4DNiftifile(fname, Fimages, structHInfo.V.mat)
fprintf('done.\n')

