%--------------------------------------------------------------------------
% Sample Matlab code for creating images with hue and alpha color-mapping.
%
%  Notes for using this code:
%  You must have OpenGL available on your system to use transparency (alpha).
%  When rendering transparency MATLAB automatically uses OpenGL if it is
%  available. If it is not available, transparency will not display.
%  See the figure property RendererMode for more information.
%
% EA Allen August 30, 2011
% eallen@mrn.org
%--------------------------------------------------------------------------
close all;clear;clc

%% 1. Load the AOD_data.mat file with sample data from the fMRI AOD experiment
%--------------------------------------------------------------------------

addpath '/Users/xli77/Dropbox (GaTech)/MISA-pytorch/figures/ISBI2022SIVA/dualcodeExample'
load '/Users/xli77/Dropbox (GaTech)/MISA/results/SIVA/fixedSubspace/mask/gmMask_TPM_thrp2_fractc8_3mm.mat'
load '/Users/xli77/Dropbox (GaTech)/MISA/results/SIVA/fixedSubspace/mask/MNI152_T1_3mm_brain.mat'

% For a single axial slice (Z = 2 mm) of data, you should have:
% spatial_map: 'Difference between Novel and Standard betas averaged over 28 subjects'
% Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
% Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
% Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
%--------------------------------------------------------------------------

%% 2. Set some defaults that will affect the appearance of the image
%--------------------------------------------------------------------------

% AY_median = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_median_ukb_3d.mat').v3d;
% label_list = {'all', 'young', 'old', 'male', 'female'};
% AY_corr = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_group_corr_ukb.mat').AY_corr;
% AY_corr_3d = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_group_corr_ukb_3d.mat').v3d;
% pvalue_3d = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_group_pvalue_ukb_3d.mat').v3d;

AY_median = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_median_sz_interaction_3d.mat').v3d;
label_list = {'all', 'young_hc', 'old_hc', 'young_pt', 'old_pt'};
% AY_median = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_median_sz_3d.mat').v3d;
% label_list = {'all', 'young', 'old', 'hc', 'pt'};
AY_corr = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_group_corr_sz_interaction.mat').AY_corr;
AY_corr_3d = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_group_corr_sz_interaction_3d.mat').v3d;
pvalue_3d = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY_group_pvalue_sz_interaction_3d.mat').v3d;

trim = 4;
num_voxels = 44318;

%%
for m = 1:2
    for s = 1:5
        for j = 2:5
            corr_v3d = abs(squeeze(pvalue_3d(s,j,:,:,:)));
            thr = 0.01/44318; %prctile(abs(AY_corr(s,j,:)), 50);
            binary_corr_v3d = zeros(size(corr_v3d));
            binary_corr_v3d(corr_v3d < thr) = 1;
            binary_corr_v3d(corr_v3d > thr) = 0;
            binary_corr_v3d(mask==0) = 0;
            filled_binary_corr_v3d = imfill(binary_corr_v3d,26,"holes");
            
            nhood = ones(2);
            se = strel(nhood);
            binary_corr_v3d_open = imopen(filled_binary_corr_v3d, se);
            smoothed_binary_corr_v3d = imclose(binary_corr_v3d_open, se);

            F = figure('Color', 'w', 'Position', [10 10 1000 150]);
            axes('Position', [0 0 1 1]);

            space = 5;
            total = 8;
            t = tiledlayout(1,total,'TileSpacing','compact','Padding','tight');

            v3d = squeeze(AY_median(m,s,j,:,:,:));
            absmax = max(abs(squeeze(v3d(:))));
            
            % if j == 1
            %     absmax = max(abs(squeeze(v3d(:))));
            % elseif (j == 2) || (j == 3)
            %     v4d = squeeze(AY_median(m,s,2:3,:,:,:));
            %     absmax = max(abs(squeeze(v4d(:))));
            % elseif (j == 4) || (j == 5)
            %     v4d = squeeze(AY_median(m,s,4:5,:,:,:));
            %     absmax = max(abs(squeeze(v4d(:))));
            % end

            % if j == 1
            %     absmax = max(abs(squeeze(v3d(:))));
            % elseif (j == 2) || (j == 5)
            %     v4d = squeeze(AY_median(m,s,[2,5],:,:,:));
            %     absmax = max(abs(squeeze(v4d(:))));
            % elseif (j == 3) || (j == 4)
            %     v4d = squeeze(AY_median(m,s,3:4,:,:,:));
            %     absmax = max(abs(squeeze(v4d(:))));
            % end

            % if j == 1
            %     absmax = max(abs(squeeze(v3d(:))));
            % else
            %     v4d = squeeze(AY_median(m,s,2:5,:,:,:));
            %     absmax = max(abs(squeeze(v4d(:))));
            % end

            % Set the Min/Max values for hue coding
            for i = 1:total
                slice = space*i+1+5;

                Underlay = squeeze(MNI152T1Template(:,:,slice)); %template
                Underlay = Underlay(:,end:-1:1)';
                spatial_map_slice = v3d(:,:,slice);
                spatial_map = spatial_map_slice(:,end:-1:1)';

                contour_map_slice = smoothed_binary_corr_v3d(:,:,slice);
                contour_map = contour_map_slice(:,end:-1:1)';

                % flip L and R
                Underlay = Underlay(:,end:-1:1);
                spatial_map = spatial_map(:,end:-1:1);
                contour_map = contour_map(:,end:-1:1);

                Underlay = Underlay(trim:end-trim,trim:end-trim);
                spatial_map = spatial_map(trim:end-trim,trim:end-trim);
                contour_map = contour_map(trim:end-trim,trim:end-trim);

                background_mask=1-logical(Underlay>0);
                background_mask_rgb = zeros(73-2*trim+1, 61-2*trim+1, 3);
                for ind = 1:3
                    background_mask_rgb(:,:,ind) = background_mask;
                end

                H_range = [-absmax absmax]; % The colormap is symmetric around zero

                % Set the Min/Max T-values for alpha coding
                A_range = [0 absmax];

                % Voxels with t-stat of 0 will be completely transparent;
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
                O_RGB = convert_to_RGB(spatial_map, CM_over, H_range);

                % WhiteBG = zeros(size(U_RGB));
                WhiteBG = U_RGB;
                WhiteBG(background_mask_rgb==1) = 1;

                % Plot the underlay
                nexttile
                %             layer1 = image(U_RGB); axis image; axis off;
                layer1 = image(WhiteBG); axis off; %axis image; axis off;
                hold on;
                % Now, add the Beta difference map as an overlay
                layer2 = image(O_RGB); %axis image

                % Use the T-statistics to create an alpha map (which must be in [0,1])

                alphamap = abs(spatial_map);
                alphamap(alphamap > A_range(2)) = A_range(2);
                alphamap(alphamap < A_range(1)) = 0;
                alphamap = alphamap/A_range(2);
                % alphamap = 1 - alphamap/A_range(2); % p-value

                % Adjust the alpha values of the overlay
                set(layer2, 'alphaData', alphamap);

                % Add some (black) contours to annotate nominal significance
                hold on;
                [C, CH] = contour(contour_map, 1, 'k');
                CH.LineWidth = 1;
                set(gca, 'visible', 'off');
                set(findall(gca, 'type', 'text'), 'visible', 'on');
                set(gca, 'XTick', [], 'YTick', []);

            end

            t.Padding = 'tight';
            t.TileSpacing = 'none';

            %%
            set(gcf,'InvertHardcopy','off');
            %         set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
            set(gcf,'Color',[1 1 1]); % RGB values [1 1 1] indicates white color

            exportgraphics(t,['/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/figures/AY_spatial_map/sz/modality',num2str(m),'_subspace',num2str(s),'_',label_list{j},'.png'],'BackgroundColor','w');
            %--------------------------------------------------------------------------

            %% 4. Create a 2D colorbar for the dual-coded overlay
            %--------------------------------------------------------------------------
            G = figure('color', 'w', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.05, 0.3]);
            x = linspace(A_range(1), A_range(2), 256);
            % x represents the range in alpha (abs(t-stats))
            y = linspace(H_range(1), H_range(2), size(CM_over,1));
            % y represents the range in hue (beta weight difference)
            [X,Y] = meshgrid(x,y); % Transform into a 2D matrix
            imagesc(x,y,Y); axis xy; % Plot the colorbar
            %         set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
            %         set(gca, 'Xcolor', 'w', 'Ycolor', 'k', 'FontSize', 17)
            set(gca, 'Xcolor', 'k', 'Ycolor', 'k', 'FontSize', 18)
            set(gca, 'YAxisLocation', 'right');
            set(gca,'yticklabel',num2str(get(gca,'ytick')','%0.3f'));
            colormap(CM_over);

            %                 hcb = colorbar;
            %                 tix = hcb.Ticks;                                            % Get Tick Values
            %                 hcb.TickLabels = compose('%0.2f',tix);                      % Set Tick Labels
            %
            alpha(X);
            alpha('scaled');
            % xlabel(alpha_label,'FontSize',12)
            % ylabel(hue_label)

            set(gcf,'InvertHardcopy','off');
            %         set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
            set(gcf,'Color',[1 1 1]); % RGB values [1 1 1] indicates white color

            saveas(gcf,['/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/figures/AY_spatial_map/sz/colorbar_modality',num2str(m),'_subspace',num2str(s),'_',label_list{j},'.png']);
        end
        close all;
    end
end

%--------------------------------------------------------------------------
