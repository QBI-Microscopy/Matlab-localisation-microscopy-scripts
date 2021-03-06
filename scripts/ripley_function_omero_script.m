 
% This script assumes that you have already segemented a region of interest
% from a localisation microscopy dataset and you have a file containing
% the XY coords from that region of interest
% DO NOT USE THIS ON A COMPLETE DATASET UNLESS YOU HAVE A SUPER COMPUTER! 
% Last update on 010515 by Daniel Matthews


clc; clear all

% check your current default group before attempting to connect
% if the file annotation you want to retrieve is not in the default
% group this script will fail

% set up
loadOmero();
W = evalin('base','whos'); %or 'base'
doesAexist = ismember('client',[W(:).name]);

if ~doesAexist
    [client,session] = connectOmero();
end

fileIds = [3639];
fileAnnotations = getFileAnnotations(session, fileIds);

num_files = numel(fileAnnotations);
downloaded = 'download_temp.csv';

% run the batch
outputdata = repmat(struct('radius',[],'ripleyK',[],'ripleyL',[]),num_files,1);
for f = 1: num_files
    
    % download the data
    getFileAnnotationContent(session,fileAnnotations(f),downloaded);
    
    % read the ROI data
    table = read_localisations(downloaded,1,3,delimiter);

    if isstruct(table)
        data = table.data;
    else
        data = table;
    end

    xcol = 1; % index of column holding x coords
    ycol = 2; % index of column holding y coords
    fcol = 3; % index of column holding first frame number

    %number of frames: set array of values for each frame
    %set limits to '[]' if not required
    %limits = [(1:4000:12001)',(4000:4000:16000)']; 
    limits = [];

    convert = false; % do you need to convert the XY coords to nm?
    camPix = 100; %nm
    coords = horzcat(data(:,xcol),data(:,ycol));
    if convert
        coords = horzcat(data(:,xcol),data(:,ycol)).*camPix;
    end

    % parameters for calculation
    % max radius for Ripley calculation in nm
    maxrad = 1000; 
    dr = 10; % step in nm
    r = 1:dr:maxrad; % distance scale array

    % bounding box
    xmin = min(coords(:,1));
    ymin = min(coords(:,2));
    xmax = max(coords(:,1));
    ymax = max(coords(:,2));
    box = [xmin, xmax, ymin, ymax];

    % the calculation
    if limits %#ok
        % a holder for the ripley data
        Kdata = zeros(length(r),length(limits));
        Ldata = zeros(length(r),length(limits));
        for frameidx = 1:length(limits);
            idx = data(:,fcol) >= limits(frameidx,1) & data(:,fcol) <= limits(frameidx,2);
            coords = coords(idx,:);
            [K,L] = ripleykfunction(coords,r,box,0);
            Kdata(:,idx) = K;
            Ldata(:,idx) = L;
        end
    else
        [Kdata,Ldata] = ripleykfunction(coords,r,box,0);
    end
    outputdata(f).radius = r;
    outputdata(f).ripleyK = Kdata;
    outputdata(f).ripleyL = Ldata;
    
    try
        delete(downloaded);
    catch
        disp('could not delete downloaded file')
    end
end
