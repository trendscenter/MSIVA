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

user = 'xli77';
addpath(['/Users/',user,'/Dropbox (GaTech)/MISA-pytorch/figures/ISBI2023/dualcodeExample']);
load(['/Users/',user,'/Dropbox (GaTech)/MISA/results/SIVA/fixedSubspace/mask/gmMask_TPM_thrp2_fractc8_3mm.mat']);
load(['/Users/',user,'/Dropbox (GaTech)/MISA/results/SIVA/fixedSubspace/mask/MNI152_T1_3mm_brain.mat']);

% For a single axial slice (Z = 2 mm) of data, you should have:
% spatial_map: 'Difference between Novel and Standard betas averaged over 28 subjects'
% Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
% Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
% Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
%--------------------------------------------------------------------------

%% 2. Set some defaults that will affect the appearance of the image
%--------------------------------------------------------------------------

% AY_median = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_median_ukb_3d.mat']).v3d;
% label_list = {'all', 'young', 'old', 'male', 'female'};
% AY_corr = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_corr_ukb.mat']).AY_corr;
% AY_corr_3d = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_corr_ukb_3d.mat']).v3d;
% pvalue_3d = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_pvalue_ukb_3d.mat']).v3d;

AY_median = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_median_sz_interaction_3d.mat']).v3d;
label_list = {'all', 'young_hc', 'old_hc', 'young_pt', 'old_pt'};
AY_corr = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_corr_sz_interaction.mat']).AY_corr;
AY_corr_3d = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_corr_sz_interaction_3d.mat']).v3d;
pvalue_3d = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_pvalue_sz_interaction_3d.mat']).v3d;

% AY_median = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_median_sz_3d.mat']).v3d;
% label_list = {'all', 'young', 'old', 'hc', 'pt'};
% AY_corr = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_corr_sz.mat']).AY_corr;
% AY_corr_3d = load(['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/mat/AY_group_corr_sz_3d.mat']).v3d;

trim = 6;
num_voxels = 44318;

%%
for m = 1:2
    for s = 1:5
        for j = 2:5
            corr = squeeze(AY_corr(s,j,:));
            prct = prctile(corr, [15, 85]);
            corr_3d = squeeze(AY_corr_3d(s,j,:,:,:));
            binary_corr_3d_pos = zeros(size(corr_3d));
            binary_corr_3d_neg = zeros(size(corr_3d));
            binary_corr_3d_10p_pos = zeros(size(corr_3d));
            binary_corr_3d_10p_neg = zeros(size(corr_3d));

            binary_corr_3d_pos(corr_3d > prct(2)) = 1;
            binary_corr_3d_neg(corr_3d < prct(1)) = 1;

            nhood = ones(2);
            se = strel(nhood);

            filled_binary_corr_3d_pos = imfill(binary_corr_3d_pos,26,"holes");
            open_binary_corr_3d_pos = imopen(filled_binary_corr_3d_pos, se);
            smoothed_binary_corr_3d_pos = imclose(open_binary_corr_3d_pos, se);
            
            filled_binary_corr_3d_neg = imfill(binary_corr_3d_neg,26,"holes");
            open_binary_corr_3d_neg = imopen(filled_binary_corr_3d_neg, se);
            smoothed_binary_corr_3d_neg = imclose(open_binary_corr_3d_neg, se);
            
            F = figure('Color', 'w', 'Position', [10 10 1200 180]);
            axes('Position', [0 0 1 1]);

            space = 5;
            total = 8;
            t = tiledlayout(1,total,'TileSpacing','compact','Padding','tight');

            v3d = squeeze(AY_median(m,s,j,:,:,:));
            % absmax = max(abs(squeeze(v3d(:))));
            
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

            if j == 1
                absmax = max(abs(squeeze(v3d(:))));
            else
                v4d = squeeze(AY_median(m,s,2:5,:,:,:));
                absmax = max(abs(squeeze(v4d(:))));
            end

            % Set the Min/Max values for hue coding
            for i = 1:total
                slice = space*i+1+5;

                Underlay = squeeze(MNI152T1Template(:,:,slice)); %template
                Underlay = Underlay(:,end:-1:1)';
                spatial_map_slice = v3d(:,:,slice);
                spatial_map = spatial_map_slice(:,end:-1:1)';

                contour_map_slice_pos = smoothed_binary_corr_3d_pos(:,:,slice);
                contour_map_pos = contour_map_slice_pos(:,end:-1:1)';

                contour_map_slice_neg = smoothed_binary_corr_3d_neg(:,:,slice);
                contour_map_neg = contour_map_slice_neg(:,end:-1:1)';

                % flip L and R
                Underlay = Underlay(:,end:-1:1);
                spatial_map = spatial_map(:,end:-1:1);
                contour_map_pos = contour_map_pos(:,end:-1:1);
                contour_map_neg = contour_map_neg(:,end:-1:1);

                Underlay = Underlay(trim:end-trim,trim:end-trim);
                spatial_map = spatial_map(trim:end-trim,trim:end-trim);
                contour_map_pos = contour_map_pos(trim:end-trim,trim:end-trim);
                contour_map_neg = contour_map_neg(trim:end-trim,trim:end-trim);

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

                % Add some contours to annotate nominal significance
                hold on;
                [Cp, CHp] = contour(contour_map_pos, 1, 'Color', 'm');
                CHp.LineWidth = 2;

                hold on;
                dark_green = [0, 0.7, 0];
                [Cn, CHn] = contour(contour_map_neg, 1, 'Color', 'k');
                CHn.LineWidth = 2;

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

            exportgraphics(t,['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/figures/AY_spatial_map/sz_r15%/modality',num2str(m),'_subspace',num2str(s),'_',label_list{j},'.png'],'BackgroundColor','w');
            %--------------------------------------------------------------------------

            %% 4. Create a 2D colorbar for the dual-coded overlay
            %--------------------------------------------------------------------------
            G = figure('color', 'w', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.05, 0.4]);
            x = linspace(A_range(1), A_range(2), 256);
            % x represents the range in alpha (abs(t-stats))
            y = linspace(H_range(1), H_range(2), size(CM_over,1));
            % y represents the range in hue (beta weight difference)
            [X,Y] = meshgrid(x,y); % Transform into a 2D matrix
            imagesc(x,y,Y); axis xy; % Plot the colorbar
            %         set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
            %         set(gca, 'Xcolor', 'w', 'Ycolor', 'k', 'FontSize', 18)
            set(gca, 'Xcolor', 'k', 'Ycolor', 'k', 'FontSize', 21);
            set(gca, 'YAxisLocation', 'right');
            set(gca,'yticklabel',num2str(get(gca,'ytick')','%0.2f'));
            colormap(CM_over);
            alpha(X);
            alpha('scaled');
            % xlabel(alpha_label,'FontSize',12)
            % ylabel(hue_label)

            set(gcf,'InvertHardcopy','off');
            %         set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
            set(gcf,'Color',[1 1 1]); % RGB values [1 1 1] indicates white color

            saveas(gcf,['/Users/',user,'/GaTech Dropbox/Xinhui Li/MSIVA/figures/IMAG2026/figures/AY_spatial_map/sz_r15%/colorbar_modality',num2str(m),'_subspace',num2str(s),'_',label_list{j},'.png']);
        end
        close all;
    end
end

%--------------------------------------------------------------------------
