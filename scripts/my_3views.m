function views = my_3views(amap3D, slices, bgv)

if ~ exist('bgv','var')
    bgv = 0;
end

sz = size(amap3D);

v1 = squeeze(amap3D(slices(1),:,:))';
v2 = squeeze(amap3D(:,slices(2),:))';
v3 = squeeze(amap3D(:,:,slices(3)))';

offset = floor((sz(2)-sz(1))/2);
offy = floor(sz(2)/4);
if bgv == 0
    views = zeros(2*offy+sum([sz(2) 2*sz(3)]),sz(2));
else
    views = bgv*ones(2*offy+sum([sz(2) 2*sz(3)]),sz(2));
end
views(offy+(1:sz(3)), :) = v1;
views(offy+sz(3)+(1:sz(3)), offset-1+(1:sz(1))) = v2;
views(offy+2*sz(3)+(1:sz(2)), offset-1+(1:sz(1))) = v3;

end