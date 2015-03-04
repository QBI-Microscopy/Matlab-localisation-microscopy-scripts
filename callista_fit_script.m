%% clear workspace and command line
clc
clear all
close all

%% load in the data

headerLines = 1; %the number of lines in the header
numCols = 13;
dataA = readLocalisations('/Users/uqdmatt2/Desktop/nrxn_01.txt',headerLines,numCols);
%dataB = readLocalisations('010.txt',headerLines,numCols);
dataB = [];
% assign the column numbers for the x and y coordinates
Xcol = 5;
Ycol = 6; %the two input files must have the same column configuration

%% Filter the localisations
%rule = 1: precision or intensity
%rule = 2: both
% rule = 1;
% limits = [10 50]; %[min max]
% dataCol = 7; % variable to be used for filtering. Use a vector for two parameters eg [1 3]
% dataA = filter_localisations(dataA,dataCol,rule,limits);
% if ~isempty(dataB)
%     dataB = filter_localisations(dataB,dataCol,rule,limits);
% end
%% set some variables:
camPixSize = 100; %pixel size on ccd in nm
originalX = 256; originalY = 256; % image size in pixels
xScale = originalX .* camPixSize;
yScale = originalY .* camPixSize;
% nmPixSizeX = xScale / image_resolution;
% nmPixSizeY = yScale / image_resolution;
%nmPixSize = sqrt(nmPixSizeX^2 + nmPixSizeY^2); % pixel size in 2D histogram
nmPixSize = 10;
image_resolution = [xScale/nmPixSize, yScale/nmPixSize]; % resolution for 2D histogram of localisation data
image_resolution(1)
%% apply channel alignment?
transformation = [];%'pre1_647_488.tform.mat'; %enter filename to apply transformation
calc_new = 0;
t_params = {transformation, calc_new};

%% set the type of correlation and the function to fit to the data
correlation = 'auto'; %'auto' for auto-correlation, 'cross' for cross-correlation
fit = 'exponential'; %name should match available fit functions
radius = 50; %in pixels

%% extract the (possibly filtered) x-y coordinates
% this step converts the coordinates pixels to nm; remove multiplication by
% camPixSize to work with data already in nm
Ax = dataA(:,Xcol);
Ay = dataA(:,Ycol);
% Axy = [Ax Ay].*camPixSize;
Axy = [Ax Ay];
if ~isempty(dataB)
    Bx = dataB(:,Xcol);
    By = dataB(:,Ycol);
    Bxy = [Bx By].*camPixSize;
else
    Bxy = [];
end

%% calculate the correlation and the fit
[correlation_data] = run_correlation_and_fit(Axy, Bxy, image_resolution, nmPixSize, t_params, [xScale yScale], correlation, fit, radius);