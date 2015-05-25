%% clear workspace and command line
clc
clear all
close all

%% load in the data

headerLines = 1; %the number of lines in the header
numCols = 14;
filename = '/Users/uqdmatt2/Desktop/automatic5.txt';
dataA = readLocalisations(filename,headerLines,numCols);

% assign the column numbers for the x and y coordinates
Xcol = 5;
Ycol = 6; %the two input files must have the same column configuration


%% set some variables:
camPixSize = 100; %pixel size on ccd in nm
originalX = 256; originalY = 256; % image size in pixels
xScale = originalX .* camPixSize;
yScale = originalY .* camPixSize;
Xrange = [0 xScale];
Yrange = [0 yScale];

nmPixSize = 50; %pixel size in 2d histogram
image_resolution = [xScale/nmPixSize, yScale/nmPixSize]; % resolution for 2D histogram of localisation data


%% extract the (possibly filtered) x-y coordinates
% this step converts the coordinates pixels to nm; remove multiplication by
% camPixSize to work with data already in nm
Ax = dataA(:,Xcol);
Ay = dataA(:,Ycol);
Axy = [Ax Ay];%.*camPixSize;


%% calculate the histogram

density(:,:,1) = hist2d(Axy ,image_resolution(1), image_resolution(2),Xrange,Yrange);
figure;hIm = imshow(density,[0 max(max(density(:,:,1)))],'XData',Xrange,'YData',Yrange); axis equal tight off;
h = imrect;
getPosition(h)

%% get the coords in the roi
xyLim = getPosition(h).*nmPixSize;
xmin = xyLim(1);
ymin = xyLim(2);
width = xyLim(3);
height = xyLim(4);
xmax = xmin + width;
ymax = ymin + height;

XPosition = Axy(:,1);
YPosition = Axy(:,2);
isROI = XPosition>xmin & XPosition<xmax & YPosition>ymin & YPosition<ymax;
X = XPosition(isROI);
Y = YPosition(isROI);

%% make new histogram of subregion
Xrange = [min(X) max(X)];
Yrange = [min(Y) max(Y)];

subregion = ceil([xmin,xmax,ymin,ymax]./nmPixSize);
roi_density = density(subregion(3):subregion(4),subregion(1):subregion(2),1);

%% create scatter histograms

sbinx = ceil((width / nmPixSize)/2); % number of bins in scatter histogram x-direction
sbiny = ceil((height / nmPixSize)/4); % number of bins in scatter histogram y-direction

figure;
subplot(3,1,1);scatter(X,Y);
set(gca,'YDir','reverse');
axis([Xrange,Yrange])
subplot(3,1,2);imagesc(roi_density,'XData',Xrange,'YData',Yrange); axis tight off;
subplot(3,1,3);hist(X,sbinx);
[xcounts,xbins] = hist(X,sbinx);
xlim([Xrange])
