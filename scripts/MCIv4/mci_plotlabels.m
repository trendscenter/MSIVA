function mci_plotlabels(F, IMinfo)
%% plotting title
for ii = 1:size(IMinfo.titlepos,1)
    Tax = axes('Parent', F, 'position', [IMinfo.titlepos(ii,1) IMinfo.titlepos(ii,2), .1 .1], 'units', 'normalized', 'color', [1 1 1]);
    axis(Tax, 'off');
    PT = text(0,0, IMinfo.datalabel, 'Parent', Tax, 'Color', [1 1 1], 'FontName', 'Arial', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

