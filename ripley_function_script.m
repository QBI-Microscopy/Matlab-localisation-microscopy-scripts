%% clear workspace and command line

 clc
 clear all
 close all

%% set some variables:

camPixSize = 100; %pixel size on ccd in nm
originalX = 256; originalY = 256; % image size in pixels
xScale = originalX .* camPixSize;
yScale = originalY .* camPixSize;
nmPixSize = 10;
image_resolutionX = xScale/nmPixSize;
image_resolutionY = yScale/nmPixSize; % resolution for 2D histogram of localisation data

radius = 1000; %in nm

limits = [(1:4000:12001)',(4000:4000:16000)']; %number of frames: set array of values for each frame
% I CHANGED 'framecol' TO 1 FOR MY TEST DATA
framecol = 9; % column containing variable to be used for filtering. 


%% load in the data

A = load('sxgfpM18DKDctlcell2.txt'); %select file
if isstruct(A)
    dataA = A.data;
else
    dataA = A;
end
%% run once to set radius
idx = dataA(:,framecol) >= limits(1,1) & dataA(:,framecol) <= limits(1,2);
data = dataA(idx,:);

%I CHANGED THE COLUMN NUMBERS HERE FROM (3,4) TO (4,5) FOR MY TEST DATA
Xcoord = data(:,3);
Ycoord = data(:,4);
Axy = horzcat(Xcoord,Ycoord).*camPixSize; %if coordinates in nm comment out camPixSize
Bxy = [];

ii = length(limits);

% calculate the correlation and the fit
frame = 0;
% [correlation_data, vq1,ripleyradius,roicoords] = run_Rip_den(Axy, Bxy, image_resolutionX, image_resolutionY, nmPixSize, [xScale yScale], radius,[],[],frame);

%fname = sprintf('clustermap0%d.tif',1); %define name here
%imwrite(uint16(vq1),fname,'tif','compression','lzw')


%% Loop over to get remaining images
% for loop over frame index 2 to length of limits: number frames - 1

for frameidx = 1:ii;
    idx = dataA(:,framecol) >= limits(frameidx,1) & dataA(:,framecol) <= limits(frameidx,2);
    data = dataA(idx,:);
    size(data)
    %I CHANGED THE COLUMN NUMBERS HERE FROM (3,4) TO (4,5) FOR MY TEST DATA
    Xcoord = data(:,3);
    Ycoord = data(:,4);
    Axy = horzcat(Xcoord,Ycoord).*camPixSize; %if coordinates in nm comment out camPixSize
    Bxy = [];
    
    density = hist2d(Axy,image_resolutionX, image_resolutionY, [0 xScale],[0 yScale]);
    fname = sprintf('2dhistogram0%d.tif',frameidx);
    imwrite(uint8(density),fname,'tif','compression','lzw')

    % calculate the correlation and the fit
%     [correlation_data, vq,~,~] = run_Rip_den(Axy, Bxy, image_resolutionX, image_resolutionY, nmPixSize, [xScale yScale], radius,ripleyradius,roicoords,frameidx);
    
%     fname = sprintf('clustermap0%d.tif',frameidx);
%     imwrite(uint16(vq),fname,'tif','compression','lzw')
    
end    
disp('finished!')
