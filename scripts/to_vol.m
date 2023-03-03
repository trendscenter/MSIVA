function out_maps_3D = to_vol(in_maps_1D, msk_3D, bgv)

if ~ exist('bgv','var')
    bgv = 0;
end

C = size(in_maps_1D,1);
sz = size(msk_3D);
if bgv == 0
    maps = zeros(prod(sz),C);
else
    maps = repmat(bgv,prod(sz),1);
end
maps(msk_3D(:),:) = in_maps_1D';
out_maps_3D = reshape(maps,[sz C]);

end