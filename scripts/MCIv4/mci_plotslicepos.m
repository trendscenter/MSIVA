function [Hf,T] = mci_plotslicepos(F, IMinfo, ind, slices, coords)

%% get the real world coordinates
RWmm{1} = coords.x(slices(1));
RWmm{2} = coords.y(slices(2));
RWmm{3} = coords.z(slices(3));

DIMletter = 'XYZ';
for ii = 1:3
    SLICEstring{ii} = sprintf('%s = %d mm', DIMletter(ii), round(RWmm{ii}));
end


%% plotting slice in mm text
for ii = 1:3 % views
    Tax = axes('Parent', F, 'position', [IMinfo.textpos(ind,ii,1) IMinfo.textpos(ind,ii,2), .1 .1], 'units', 'normalized', 'color', [1 1 1]);
    axis(Tax, 'off');
    T(ii) = text(0,0, SLICEstring{ii}, 'Parent', Tax, 'Color', [1 1 1], 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'center','VerticalAlignment', 'top');
end
