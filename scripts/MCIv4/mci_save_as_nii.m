function V = mci_save_as_nii(data, sesInfo, outpath, outname)

B = zeros(sesInfo.HInfo.DIM);
B(sesInfo.mask_ind) = data;
tempV = sesInfo.HInfo.V;
tempV.fname = fullfile(outpath,outname);
tempV = rmfield(tempV, 'pinfo');
V = spm_write_vol(tempV,B);