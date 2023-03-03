%--------------------------------------------------------------------------
% Sample Matlab code for creating images with hue and alpha color-mapping.
%
%  Notes for using this code:
%  You must have OpenGL available on your system to use transparency (alpha).
%  When rendering transparency MATLAB automatically uses OpenGL if it is
%  available. If it is not available, transparency will not display.
%  See the figure property RendererMode for more information.
%
% Original code was written by EA Allen August 30, 2011 eallen@mrn.org
% Current code was written by Xinhui Li March 1, 2023 xli993@gatech.edu
%--------------------------------------------------------------------------

%% 1. Load data
%--------------------------------------------------------------------------
close all;clear;clc

load /Users/xli77/Documents/MSIVA/MISA-data/dualmap/CCA_vol3d_10scv.mat % UKB spatial map
% load /Users/xli77/Documents/MSIVA/MISA-data/dualmap/CCA_vol3d_sz_10scv.mat % patient spatial map
load /Users/xli77/Documents/MSIVA/MISA-data/dualmap/gmMask_TPM_thrp2_fractc8_3mm.mat % GM mask
load /Users/xli77/Documents/MSIVA/MISA-data/dualmap/MNI152_T1_3mm_brain.mat % T1 background

%--------------------------------------------------------------------------

%% 2. Set some defaults that will affect the appearance of the image
%--------------------------------------------------------------------------

for scv = 1:10
    for m = 1:2

        F = figure('Color', 'w', 'Position', [10 10 1000 150]);
        axes('Position', [0 0 1 1]);

        space = 5;
        total = 8;
        t = tiledlayout(1,total,'TileSpacing','none','Padding','tight');

        if mod(scv,2)==0
            v = vol3d(m,scv-1:scv,:,:,:);
        else
            v = vol3d(m,scv:scv+1,:,:,:);
        end

        absmax = max(abs(squeeze(v(:))));

        % Set the Min/Max values for hue coding
        for i = 1:total

            slice = space*i+1+5;

            background_mask_rgbap_N_S = squeeze(vol3d(m,scv,:,:,slice));

            Underlay = squeeze(MNI152T1Template(:,:,slice));

            background_mask_rgbap = background_mask_rgbap_N_S(mask(:,:,slice)>0);
            Q75 = quantile(background_mask_rgbap,0.75);
            Q25 = quantile(background_mask_rgbap,0.25);

            Pmap_N_S = zeros(size(background_mask_rgbap_N_S));
            Pmap_N_S(abs(background_mask_rgbap_N_S)>max(abs(Q25),abs(Q75))) = 1; %max(Q25,Q75)

            Underlay = Underlay(:,end:-1:1)';
            background_mask_rgbap_N_S = background_mask_rgbap_N_S(:,end:-1:1)';
            Pmap_N_S = Pmap_N_S(:,end:-1:1)';

            % flip L and R
            Underlay = Underlay(:,end:-1:1);
            background_mask_rgbap_N_S = background_mask_rgbap_N_S(:,end:-1:1);
            Pmap_N_S = Pmap_N_S(:,end:-1:1);

            background_mask=1-logical(Underlay>0);
            background_mask_rgb = zeros(73, 61, 3); 
            for ind = 1:3
                background_mask_rgb(:,:,ind) = background_mask;
            end

            H_range = [-absmax absmax]; % The colormap is symmetric around zero

            % Set the Min/Max T-values for alpha coding
            A_range = [0 absmax];

            % Voxels with t-stats of 0 will be completely transparent;
            % voxels with t-stat magnitudes greater or equal than 5 will be opaque.

            % Set the labels for the colorbar
            hue_label = '';
            alpha_label = '';

            % Choose a colormap for the underlay
            CM_under = gray(256);

            % Choose a colormap for the overlay
            CM_over = jet(256);
            %--------------------------------------------------------------------------

            %% 3. Do the actual plotting
            %--------------------------------------------------------------------------
            % Make a figure and set of axes

            % Transform the underlay and beta map to RGB values, based on specified colormaps
            % See function convert_to_RGB() for more information
            U_RGB = convert_to_RGB(Underlay, CM_under);
            O_RGB = convert_to_RGB(background_mask_rgbap_N_S, CM_over, H_range);

            WhiteBG = U_RGB;
            WhiteBG(background_mask_rgb==1) = 1;

            % Plot the underlay
            nexttile
            layer1 = image(WhiteBG); axis off; %axis image; axis off;
            hold on;
            % Now, add the Beta difference map as an overlay
            layer2 = image(O_RGB); %axis image

            % Use the T-statistics to create an alpha map (which must be in [0,1])
            alphamap = abs(background_mask_rgbap_N_S);
            alphamap(alphamap > A_range(2)) = A_range(2);
            alphamap(alphamap < A_range(1)) = 0;
            alphamap = alphamap/A_range(2);

            % Adjust the alpha values of the overlay
            set(layer2, 'alphaData', alphamap);

            % Add some (black) contours to annotate nominal significance
            hold on;
            [C, CH] = contour(Pmap_N_S, 1, 'k');
            set(gca, 'visible', 'off');
            set(findall(gca, 'type', 'text'), 'visible', 'on');
            set(gca, 'XTick', [], 'YTick', []);

        end

        t.Padding = 'tight';
        t.TileSpacing = 'none';

        %%
        set(gcf,'InvertHardcopy','off');
        set(gcf,'Color',[1 1 1]); % RGB values [1 1 1] indicates white color

%         exportgraphics(t,['/Users/xli77/Documents/MSIVA/figures/ISBI2023SIVA/Source',num2str(scv),'M',num2str(m),'.png'],'BackgroundColor','w');

        %--------------------------------------------------------------------------

        %% 4. Create a 2D colorbar for the dual-coded overlay
        %--------------------------------------------------------------------------
        G = figure('color', 'w', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.04, 0.35]);
        x = linspace(A_range(1), A_range(2), 256);
        % x represents the range in alpha (abs(t-stats))
        y = linspace(H_range(1), H_range(2), size(CM_over,1));
        % y represents the range in hue (beta weight difference)
        [X,Y] = meshgrid(x,y); % Transform into a 2D matrix
        imagesc(x,y,Y); axis xy; % Plot the colorbar
        set(gca, 'Xcolor', 'w', 'Ycolor', 'k', 'FontSize', 17)
        set(gca, 'YAxisLocation', 'right');
        set(gca,'yticklabel',num2str(get(gca,'ytick')','%0.2f'));
        colormap(CM_over);

        set(gcf,'InvertHardcopy','off');
        set(gcf,'Color',[1 1 1]); % RGB values [1 1 1] indicates white color

%         saveas(gcf, ['/Users/xli77/Documents/MSIVA/figures/ISBI2023SIVA/colorbar_Source',num2str(scv),'M',num2str(m),'.png']); % save as .png file

    end
end
%--------------------------------------------------------------------------
