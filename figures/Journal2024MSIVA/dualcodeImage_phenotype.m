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
ver = 4;
n_var = 25;
addpath(['/Users/',user,'/Dropbox (GaTech)/MISA-pytorch/figures/ISBI2022SIVA/dualcodeExample']);
gm_path = ['/Users/',user,'/Dropbox (GaTech)/MISA/results/SIVA/fixedSubspace/mask/gmMask_TPM_thrp2_fractc8_3mm.mat'];
brain_path = ['/Users/',user,'/Dropbox (GaTech)/MISA/results/SIVA/fixedSubspace/mask/MNI152_T1_3mm_brain.mat'];
load(gm_path);
load(brain_path);
dir = ['/Users/',user,'/Dropbox (GaTech)/MISA-pytorch/figures/SIVA/v',num2str(ver),'/phenotype_map_',num2str(n_var),'var/'];
p_path = [dir,'p_corrected3d.mat'];
reject_path = [dir,'rejected3d.mat'];
load(p_path);
load(reject_path);

coeff_var = load('/Users/xli77/GaTech Dropbox/Xinhui Li/MSIVA/figures/Journal2024MSIVA/mat/AY5_coeff_var_ukb_3d.mat').v3d;
absmax_alpha = max(coeff_var(:));

% For a single axial slice (Z = 2 mm) of data, you should have:
% spatial_map: 'Difference between Novel and Standard betas averaged over 28 subjects'
% Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
% Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
% Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
%--------------------------------------------------------------------------

%% 2. Set some defaults that will affect the appearance of the image
%--------------------------------------------------------------------------
absmax = 0;
for pheno = 1:n_var
    for ss = 1:14
        r3d = squeeze(load([dir,'phenotype_map',num2str(pheno),'_r3d.mat']).r3d(ss,:,:,:));
        tmp = max(abs(r3d(:)));
        if tmp > absmax
            absmax = tmp;
        end
    end
end

%%
for pheno = [4,13,18]
    for ss = 5
        F = figure('Color', 'w', 'Position', [10 10 1000 150]);
        axes('Position', [0 0 1 1]);

        space = 5;
        total = 8;
        t = tiledlayout(1,total,'TileSpacing','none','Padding','tight');

        r3d = squeeze(load([dir,'phenotype_map',num2str(pheno),'_r3d.mat']).r3d(ss,:,:,:));
        p3d = squeeze(load([dir,'phenotype_map',num2str(pheno),'_p3d.mat']).p3d(ss,:,:,:));
        % absmax = max(abs(r3d(:)));
        
        rejected = rejected3d(:,:,:,ss,pheno);
        
        % Set the Min/Max values for hue coding
        for i = 1:total
            slice = space*i+1+5;

            Underlay = squeeze(MNI152T1Template(:,:,slice)); %template
            Underlay = Underlay(:,end:-1:1)';
            v = r3d(:,:,slice);
            spatial_map = v(:,end:-1:1)';
            
            Pmap_N_S = rejected(:,:,slice);
            Pmap_N_S(mask(:,:,slice)==0) = 0;
            Pmap_N_S = Pmap_N_S(:,end:-1:1)';

            coeff_var_map = coeff_var(:,:,slice);
            coeff_var_map = coeff_var_map(:,end:-1:1)';

            % flip L and R
            Underlay = Underlay(:,end:-1:1);
            spatial_map = spatial_map(:,end:-1:1);
            Pmap_N_S = Pmap_N_S(:,end:-1:1);
            coeff_var_map = coeff_var_map(:,end:-1:1);
            
            background_mask=1-logical(Underlay>0);
            background_mask_rgb = zeros(73, 61, 3);
            for ind = 1:3
                background_mask_rgb(:,:,ind) = background_mask;
            end

            H_range = [-absmax absmax]; % The colormap is symmetric around zero

            % Set the Min/Max T-values for alpha coding
            A_range = [0 absmax_alpha];

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
            % alphamap = ones(size(spatial_map));
            alphamap = abs(coeff_var_map);
            alphamap(alphamap > A_range(2)) = A_range(2);
            alphamap(alphamap < A_range(1)) = 0;
            alphamap = alphamap/A_range(2);

            % Adjust the alpha values of the overlay
            set(layer2, 'alphaData', alphamap);

            % Add some (black) contours to annotate nominal significance
            %             hold on;
            [C, CH] = contour(Pmap_N_S, 1, 'k');
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

        exportgraphics(t,['/Users/',user,'/Dropbox (GaTech)/MISA-pytorch/figures/SIVA/v',num2str(ver),'/phenotype_map_',num2str(n_var),'var_fdr_corrected/phenotype',num2str(pheno),'_predictor',num2str(ss),'.png'],'BackgroundColor','w')
        %--------------------------------------------------------------------------

        %% 4. Create a 2D colorbar for the dual-coded overlay
        %--------------------------------------------------------------------------
        %         G = figure('color', 'k', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.06, 0.35]);
        G = figure('color', 'w', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.059, 0.35]);
        x = linspace(A_range(1), A_range(2), 256);
        % x represents the range in alpha (abs(t-stats))
        y = linspace(H_range(1), H_range(2), size(CM_over,1));
        % y represents the range in hue (beta weight difference)
        [X,Y] = meshgrid(x,y); % Transform into a 2D matrix
        imagesc(x,y,Y); axis xy; % Plot the colorbar
        set(gca, 'Xcolor', 'k', 'Ycolor', 'k', 'FontSize', 17);
        set(gca, 'YAxisLocation', 'right');
        set(gca,'yticklabel',num2str(get(gca,'ytick')','%0.2f'));
        colormap(CM_over);

        alpha(X);
        alpha('scaled');
        % xlabel(alpha_label,'FontSize',14)
        % ylabel(hue_label)

        set(gcf,'InvertHardcopy','off');
        %         set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
        set(gcf,'Color',[1 1 1]); % RGB values [1 1 1] indicates white color

        saveas(gcf,['/Users/',user,'/Dropbox (GaTech)/MISA-pytorch/figures/SIVA/v',num2str(ver),'/phenotype_map_',num2str(n_var),'var_fdr_corrected/colorbar_phenotype',num2str(pheno),'_predictor',num2str(ss),'.png'])
    end
    % close all;
end

%--------------------------------------------------------------------------
