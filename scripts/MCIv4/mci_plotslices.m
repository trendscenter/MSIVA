function Hf = mci_plotslices(F, IMinfo, ind, data, slices, cm, coords)


%% plotting composite maps

Hf(1) = subplot('position', squeeze(IMinfo.sgrid(ind,1,:))');
hold off
imagesc(squeeze(data(slices(1),:,:))'); axis image;  set(gca, 'clim', [1, length(cm)]); axis xy; axis off

Hf(2) = subplot('position', squeeze(IMinfo.sgrid(ind,2,:))');
hold off
imagesc(squeeze(data(:,slices(2),:))'); axis image;  set(gca, 'clim', [1, length(cm)]); axis xy; axis off

Hf(3) = subplot('position', squeeze(IMinfo.sgrid(ind,3,:))');
hold off
imagesc(squeeze(data(:,:,slices(3)))'); axis image;  set(gca, 'clim', [1, length(cm)]); axis xy; axis off
colormap(cm);






