%% clear workspace and command line

 clc
 clear all
 close all

%% set some variables:

camPixSize = 95; %pixel size on ccd in nm
originalX = 180; originalY = 180; % image size in pixels
xScale = originalX .* camPixSize;
yScale = originalY .* camPixSize;
nmPixSize = 10;
image_resolutionX = xScale/nmPixSize;
image_resolutionY = yScale/nmPixSize; % resolution for 2D histogram of localisation data

%camPixSize = 100; %pixel size on ccd in nm
%originalX = 326; originalY = 282; % image size in pixels
%image_resolutionX = 3260;   % resolution for 2D histogram of localisation data
%image_resolutionY = 2820;   % resolution for 2D histogram of localisation data
%nmPixSizeX = xScale / image_resolutionX;
%nmPixSizeY = yScale / image_resolutionY;
%nmPixSize = sqrt(nmPixSizeX^2 + nmPixSizeY^2); 

radius = 50; %in pixels

limits = [(0:250:1000)',(2000:250:3000)']; %number of frames: set array of values for each frame
framecol = 1; % column containing variable to be used for filtering. 


%% load in the data

A = load('test_data/488_011.txt'); %select file
if isstruct(A)
    dataA = A.data;
else
    dataA = A;
end

%% run once to set radius
idx = dataA(:,framecol) >= limits(1,1) & dataA(:,framecol) <= limits(1,2);
data = dataA(idx,:);

Xcoord = data(:,4);
Ycoord = data(:,5);
Axy = horzcat(Xcoord,Ycoord).*camPixSize; %if coordinates in nm comment out camPixSize
Bxy = [];

ii = length(limits);

% calculate the correlation and the fit
frame = 1;
[correlation_data, vq1,ripleyradius,roicoords] = run_Rip_den(Axy, Bxy, image_resolutionX, image_resolutionY, nmPixSize, [xScale yScale], radius,[],[],frame);

fname = sprintf('clustermap0%d.tif',1); %define name here
imwrite(uint16(vq1),fname,'tif','compression','lzw')


%% Loop over to get remaining images
% for loop over frame index 2 to length of limits: number frames - 1

for frameidx = 2:ii;
    idx = dataA(:,framecol) >= limits(frameidx,1) & dataA(:,framecol) <= limits(frameidx,2);
    data = dataA(idx,:);

    Xcoord = data(:,4);
    Ycoord = data(:,5);
    Axy = horzcat(Xcoord,Ycoord).*camPixSize; %if coordinates in nm comment out camPixSize
    Bxy = [];

    % calculate the correlation and the fit
    [correlation_data, vq,~,~] = run_Rip_den(Axy, Bxy, image_resolutionX, image_resolutionY, nmPixSize, [xScale yScale], radius,ripleyradius,roicoords,frameidx);
    
    fname = sprintf('clustermap0%d.tif',frameidx);
    imwrite(uint16(vq),fname,'tif','compression','lzw')
    
end    
disp('finished!')
