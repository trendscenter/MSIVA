function [sgrid, textpos, titlepos, cbarpos] = mci_makesubplotgrid_cluster(ncomps, colorbarmargin, vdim, slicedir)

%% vdim = [xdim, ydim, zdim];
%% sgrid is ncomp x dim x [left bottom width height]
% colorbarmargin = .1;
% vdim = size(MAP{1});
% slicedir = 'h' (horizontal subplots) or 'v' (vertical subplots)

plotfrac = 1;

fullw = 1-colorbarmargin;

if strcmp(slicedir, 'h')
    hpixels = vdim(2) + 2*vdim(1);
else
    hpixels = ncomps*max(vdim);
end
hplot1 = (vdim(2)/hpixels)*(plotfrac)*fullw;
hplot2 = (vdim(1)/hpixels)*(plotfrac)*fullw;
hplot3 = (vdim(1)/hpixels)*(plotfrac)*fullw;

hspace = (fullw-(hplot1+hplot2+hplot3))/3;


hspace1 = (fullw-ncomps*hplot1)/(ncomps+1);
hspace2 = (fullw-ncomps*hplot2)/(ncomps+1);
hspace3 = (fullw-ncomps*hplot3)/(ncomps+1);

if strcmp(slicedir, 'h')
    vpixels = ncomps*max(vdim);
else
    vpixels = 2*vdim(3)+vdim(2);
end

vplot1 = (vdim(3)/vpixels)*(plotfrac)*fullw;
vplot2 = (vdim(3)/vpixels)*(plotfrac)*fullw;
vplot3 = (vdim(2)/vpixels)*(plotfrac)*fullw;

vspace = (1-(vplot1+vplot2+vplot3))/3;

vspace1 = (1-ncomps*vplot1)/(ncomps+1);
vspace2 = (1-ncomps*vplot2)/(ncomps+1);
vspace3 = (1-ncomps*vplot3)/(ncomps+1);

%% SUBPLOT('position',[left bottom width height])
if strcmp(slicedir, 'h')
    for ii =1:ncomps
        sgrid(ii,1,:) = [hspace,  (ncomps-ii+1)*(vspace1)+(ncomps-ii)*vplot1, hplot1, vplot1];
        sgrid(ii,2,:) = [2*hspace+hplot1, (ncomps-ii+1)*(vspace2)+(ncomps-ii)*vplot2, hplot2, vplot2];
        sgrid(ii,3,:) = [3*hspace+hplot1+hplot2, (ncomps-ii+1)*(vspace3)+(ncomps-ii)*vplot3, hplot3, vplot3];
    end
else
    for ii =1:ncomps
        sgrid(ii,1,:) = [(ii)*hspace1+(ii-1)*hplot1, 3*vspace+vplot3+vplot2,    hplot1, vplot1];
        sgrid(ii,2,:) = [(ii)*hspace2+(ii-1)*hplot2, 2*vspace+vplot3,           hplot2, vplot2];
        sgrid(ii,3,:) = [(ii)*hspace3+(ii-1)*hplot3, vspace,                    hplot3, vplot3];
    end
end


%% text info
%textpos is [x x x; y y y]
%ty = (sgrid(end-1,1,2) + (sgrid(end,1,2)+sgrid(end,1,4)))/2;
%% slice info 

avgvspace = ((vspace1+vspace2+vspace3)/3);
avgvplot = (vplot1+vplot2+vplot3)/3;

textpos = zeros(ncomps,3,2);
for ii=1:ncomps
    textpos(ii,1,1) = sgrid(ii,1,1)+0.5*sgrid(ii,1,3); %xpos1
    textpos(ii,2,1) = sgrid(ii,2,1)+0.5*sgrid(ii,2,3); %xpos2
    textpos(ii,3,1) = sgrid(ii,3,1)+0.5*sgrid(ii,3,3); %xpos3
    if strcmp(slicedir, 'h')
        textpos(ii,:,2) = sgrid(ii,3,2)+avgvspace*0.3;   %ypos
    else
        
        textpos(ii,1,2) = sgrid(ii,1,2)+avgvspace*0.1;   %ypos1
        textpos(ii,2,2) = sgrid(ii,2,2)+avgvspace*0.1;   %ypos2
        textpos(ii,3,2) = sgrid(ii,3,2)+avgvspace*0.1;   %ypos3
    end
end


% ty = avgvspace/2;
% txshift = .5*[hplot1 hplot2 hplot3];
% textpos = [txshift+[hspace, 2*hspace+hplot1,  3*hspace+hplot1+hplot2]; ty ty ty];

% %title pos
titlepos = zeros(1,2); % for just the first cluster
for ii=1
    titlepos(ii,1) = sgrid(1,2,1)+hplot2/2; %fullw/2; %xpos
    titlepos(ii,2) = avgvspace*(ncomps-ii+1) + (ncomps-ii+1)*(avgvplot)-avgvspace*0.5;    %ypos
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% colobar bar info
for ii = 1 % only print 1 colorbar
    xm = colorbarmargin/4;
    height = vplot3*.8;
    bottom = sgrid(ii,3,2)+avgvspace*0.3;
    cbarpos{ii} = [fullw, bottom,  xm, height];
end


