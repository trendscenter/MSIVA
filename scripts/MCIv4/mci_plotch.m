function mci_plotch(F, Hf, IMinfo, data, slices)

%% draw cross hairs
if ~isempty(IMinfo.crosshair_color)
    for ii = 1:3
        axes(Hf(ii));
        hold on
        ch_margin = 0.05;
        yline = [slices(2), slices(2)];
        xline = [slices(1), slices(1)];
        zline = [slices(3),slices(3)];
        xlim =  [1+size(data,1)*ch_margin, size(data,1)*(1-ch_margin)];
        ylim =  [1+size(data,2)*ch_margin, size(data,2)*(1-ch_margin)];
        zlim =  [1+size(data,3)*ch_margin, size(data,3)*(1-ch_margin)];

        if ii == 1
            %vline
            L(1) = plot(yline, zlim);
            %hline
            L(2) = plot(ylim, zline);
        elseif ii == 2
            %vline
            L(1) = plot(xline, zlim);
            %hline
            L(2) = plot(xlim, zline);
        elseif ii == 3
            %vline
            L(1) = plot(xline, ylim);
            %hline
            L(2) = plot(xlim, yline);
        end
        set(L, 'Color', IMinfo.crosshair_color);
    end
end