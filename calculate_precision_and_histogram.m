%% clear workspace and command line
clc
clear all
close all

%% load in the data

headerLines = []; %the number of lines in the header
numCols = 12;
dataA = readLocalisations('ap180_glrt35.txt',numCols,headerLines);
dataB = readLocalisations('gm1_glrt35.txt',numCols,headerLines);%readLocalizer('019 combined emitters.txt');

%% assign the column numbers for the x and y coordinates
camPixSize = 99; %pixel size on ccd in nm
NphotonsCol = 2;
SDCol = 3;
Xcol = 4;
Ycol = 5; %the two input files must have the same column configuration
BGCol = 6;
% calculate the localisation precision
loc_precisionA = ones(size(dataA,1),1).*NaN;
gainFactor = 54; %this is from France metadata
for ii = 1: size(dataA,1)
    NumPhotons = dataA(ii,NphotonsCol) / gainFactor;
    NumBGPhotons = dataA(ii,BGCol)  / gainFactor;
    loc_precisionA(ii) = precision(dataA(ii,SDCol),camPixSize,NumPhotons,NumBGPhotons);
end
dataA(:,end+1) = loc_precisionA;
if ~isempty(dataB)
    loc_precisionB = ones(size(dataB,1),1).*NaN;
    for ii = 1: size(dataB,1)
        NumPhotons = dataB(ii,NphotonsCol) / gainFactor;
        NumBGPhotons = dataB(ii,BGCol)  / gainFactor;
        loc_precisionB(ii) = precision(dataB(ii,SDCol),camPixSize,NumPhotons,NumBGPhotons);
    end
    dataB(:,end+1) = loc_precisionB;
end
%% Plot histogram
dataCol = 13; %which column should be histogrammed
bins = 100; %how many bins in the histogram
%x = linspace(0,100,bins);
variableValuesA = dataA(:,dataCol);
figure;hist(variableValuesA,bins);
if ~isempty(dataB)
    variableValuesB = dataB(:,dataCol);
    figure;hist(variableValuesB,bins);
end
