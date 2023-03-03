function [FH, S] = mci_plotcomps(MAP, structuralImage, maxncluster, minclustersize, threshold, Gmax, absflag, coords, datalabel, cbarlabel, crosshair_color, CM, slices, slicedir)


%% TO make LEFT appear on the LEFT
MAP = flipdim(MAP,1);

%% threshold
if threshold > 0
    if absflag
        MAP(find(abs(MAP) < threshold)) = 0;
    else
        MAP(find(MAP < threshold)) = 0;
    end
end

%% Mask for which functional data should be plotted
MASK = find(MAP ~= 0);

if ischar(maxncluster)
    interactive = 1;
    maxncluster = 1;
else
    interactive = 0;
end

%% where is the middle?
midshift = 5; %number of slices to move over to not plot exactly on midline 
%(may need to be adjusted based on structural dimensions)
[mv, midslice] = (min(abs(coords.x-0)));

%% if no slices are provided, find them
if ischar(slices)
    if maxncluster == 1
        B = zeros(size(MAP));
        B(MASK) = MAP(MASK);
        [S.x, S.y, S.z] = get_slices(B, slices);
        if abs(S.x - midslice) < midshift
            S.x = midslice + midshift;
        end
        ncluster = 1;
    else
        %% Clustering
        % binarize the thresholded map
        B = zeros(size(MAP));
        B(MASK) = 1;
        fprintf('clustering...\n')
        % use bwlabel to find clusters exceeding size
        D = mci_getclusters(B, minclustersize);
        fprintf('done\n')

        if isempty(D)
            disp('No clusters exceeding minumum cluster size')
            FH = [];
            S = [];
            return
        end
        %% sort clusters by size
        [ D.n, sIND] =sort(D.n, 'descend');
        D.mask = D.mask(sIND);

        %% number of clusters to plot
        ncluster = min([length(D.n), maxncluster]);

        %% Get the coordinates for each cluster
        for cc = 1:ncluster
            B = zeros(size(MAP));
            B(D.mask{cc}) = MAP(D.mask{cc});
            [S.x(cc), S.y(cc), S.z(cc)] = get_slices(B, slices);
            if abs(S.x(cc) - midslice) < midshift
                S.x(cc) = midslice + midshift;
            end
        end
    end % of 1 cluster
else % slices are provided
    S = slices;
    ncluster = length(S.x);
end


%% Limits of functional/structural data
Gmin = threshold;
Grange = Gmax-Gmin;
sLIM = [min(structuralImage(:)), max(structuralImage(:))];
if absflag
    Fun_range = [-Gmax Gmax];
else
    Fun_range = [0 Gmax];
end

cm = makecolormap(CM, Fun_range, threshold, absflag);
colorlength = length(cm);
cm = [cm;  gray(colorlength);];

FH = figure; set(FH, 'Color', 'k');
scrsz = get(0,'ScreenSize');
if strcmp(slicedir, 'h')
    height = min(scrsz(4)*0.95, maxncluster*.26*scrsz(4));
    bottom = scrsz(4)/2 - height/2;
    set(FH, 'OuterPosition',[scrsz(3)/4, bottom, scrsz(3)/3, height]);
else
    height = scrsz(4)*2/3;
    bottom = scrsz(4)/2 - height/2;
    width =  min(scrsz(3)*0.26, maxncluster*.15*scrsz(3));
    set(FH, 'OuterPosition',[scrsz(3)/4, bottom, width, height]);
end
%[left bottom width height]

colorbarmargin = 0.075;
[IMinfo.sgrid, IMinfo.textpos, IMinfo.titlepos, IMinfo.cbarpos] = mci_makesubplotgrid_cluster(maxncluster, colorbarmargin, size(MAP), slicedir);
IMinfo.minval = Fun_range(1);
IMinfo.maxval = Fun_range(2);
IMinfo.datalabel = datalabel;
IMinfo.cbarlabel = cbarlabel;
IMinfo.crosshair_color = crosshair_color;

%% Remap to colorscale
[sfdata, ssdata] = remap_to_colorscale(MAP, structuralImage, Fun_range, sLIM, colorlength);
pdata = ssdata;
%% Fill in functional on structural
pdata(MASK) = sfdata(MASK);
%% PLOT!
for cc = 1:ncluster    
    axH = mci_plotslices(FH, IMinfo, cc, pdata, [S.x(cc), S.y(cc), S.z(cc)], cm, coords);  
    mci_plotch(FH, axH, IMinfo, pdata, [S.x(cc), S.y(cc), S.z(cc)]);
    mci_plotslicepos(FH, IMinfo, cc, [S.x(cc), S.y(cc), S.z(cc)], coords);
    mci_plotlabels(FH, IMinfo);
end
mci_plotcolorbar(FH, IMinfo, 1, cm);
%% Interactive viewing!

if interactive

    fprintf('\tLeft-click to update slices.  Right-click when done.\n')
    [thispoint(1), thispoint(2), pickslice] = ginput(1);
    while pickslice == 1
        thisAx = gca;
        thisAx = find(thisAx == axH);

        if thisAx <= 3
            if thisAx == 1
                S.x = S.x;
                S.y = round(thispoint(1));
                S.z = round(thispoint(2));
            elseif thisAx ==2
                S.y = S.y;
                S.x = round(thispoint(1));
                S.z = round(thispoint(2));
            else
                S.z = S.z;
                S.x = round(thispoint(1));
                S.y = round(thispoint(2));
            end
            %update slices
            axH = mci_plotslices(FH, IMinfo, cc, pdata, [S.x, S.y, S.z], cm, coords);
            mci_plotslicepos(FH, IMinfo, cc, [S.x(cc), S.y(cc), S.z(cc)], coords);
        else
            fprintf('Please click on sagittal, coronal, or axial views to update slices, or hit enter to exit.\n')
        end
       
        [thispoint(1), thispoint(2), pickslice] = ginput(1);

    end
    axH = mci_plotslices(FH, IMinfo, cc, pdata, [S.x(cc), S.y(cc), S.z(cc)], cm, coords);
    mci_plotch(FH, axH, IMinfo, pdata, [S.x(cc), S.y(cc), S.z(cc)]);
    mci_plotslicepos(FH, IMinfo, cc, [S.x(cc), S.y(cc), S.z(cc)], coords);
    mci_plotlabels(FH, IMinfo);
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%sub-functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, z] = get_slices(B, slice_method)

switch lower(slice_method)
    case {'dvar'}
        temp = squeeze(var(B,[],3));
        [junk, x_sorted] = sort(var(temp,[],2), 'descend');
        [junk, y_sorted] = sort(var(temp,[],1),'descend' );
        temp = squeeze(var(B,[],2));
        [junk, z_sorted] = sort(var(temp,[],1), 'descend');
        x = x_sorted(1);
        y = y_sorted(1);
        z = z_sorted(1);
    case {'dmax'}
        [maxval, maxind] = max(abs(B(:)));
        [x, y, z] = ind2sub(size(B),maxind);
    case {'dsum'}
        temp = squeeze(sum(B,3));
        [junk, x_sorted] = sort(sum(temp,2), 'descend');
        [junk, y_sorted] = sort(sum(temp,1),'descend' );
        temp = squeeze(sum(B,2));
        [junk, z_sorted] = sort(sum(temp,1), 'descend');
        x = x_sorted(1);
        y = y_sorted(1);
        z = z_sorted(1);
end


function CM = makecolormap(cm, range, threshold, absflag)

%maxlen = 64;
nskip = 1;%ceil(length(cm)/maxlen);
% uncomment lines below to pad the colormap with zeros, scaled with the
% threshold; otherwise, colormap begins with the lowest value passing the
% threshold
if absflag
    cm = cm(1:nskip:end,:);
    cm1 = cm(1:floor(end/2),:);
    cm2 = cm(floor(end/2)+1:end,:);
    %zerolength = 1;
    zerolength = round(threshold*length(cm1)/(range(2)-threshold));
    CM = [cm1; zeros(zerolength,3); cm2];
else
    cm = cm(1:nskip:end,:);
    %zerolength = 1;
    zerolength = floor(threshold*length(cm)/(range(2)-threshold));
    CM = [zeros(zerolength,3); cm];
end

function [sfdata, ssdata] = remap_to_colorscale(fdata, sdata, fLIM, sLIM, colorlength)
sfdata = scale_data(fdata, fLIM, [1, colorlength]);
ssdata = scale_data(sdata, sLIM, [colorlength+1, 2*colorlength]);

function sd = scale_data(d, oldLIM, newLIM)
old_range = diff(oldLIM);
new_range = diff(newLIM);
sd = (((d-oldLIM(1))./old_range)./(1/new_range)) + newLIM(1);

