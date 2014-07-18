function [corrData,vq,selectedRadius,roicoords] = run_Rip_den(varargin) %(Axy, Bxy, image_resolutionX, image_resolutionY, nmPixSize, [xScale yScale], correlation, fit, radius);

    channel1 = varargin{1};
    channel2 = varargin{2};
    resX = varargin{3};
    resY = varargin{4};
    nmpixSize = varargin{5};
    range = varargin{6};
    Xrange = [0 range(1)];
    Yrange = [0 range(2)];
    maxrad1 = varargin{7};
    ripleyradius = varargin{8};
    roicoords = varargin{9};
    frame = varargin{10};
    
    N = {};
    data = {};
    N{1} = size(channel1,1);
    if isempty(channel2)
        numChannels = 1;  
        data{1} = channel1;
        data{2} = [];
    else
        numChannels = 2;
        N{2} = size(channel2,1);
        if ~isempty(t_params{1}) && t_params{2} == 0
            tformData = load(t_params{1});
            TFORM = tformData.TFORM;
            data{1} = tformfwd(TFORM,channel1);
        elseif isempty(t_params{1}) && t_params{2} == 0
            data{1} = channel1;
        elseif isempty(t_params{1}) && t_params{2} == 1
            [in_points,base_points,~,TFORM] = transformChannels(nmpixSize,density);
            data{1} = tformfwd(TFORM,channel1);
        end
        data{2} = channel2;
    end
    
    if frame == 1
%        %calculate histograms
%         if ~isempty(channel2)
%             density(:,:,1) = hist2d(channel1,resX, resY,Xrange,Yrange);
%             density(:,:,2) = hist2d(channel2,resX, resY,Xrange,Yrange);
%             density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
%         else
%             density = hist2d(channel1,resX, resY,Xrange,Yrange);
%         end    

            %calculate histograms
        if numChannels == 2
            density(:,:,1) = hist2d(data{1},resX, resY,Xrange,Yrange);
            density(:,:,2) = hist2d(data{2},resX, resY,Xrange,Yrange);
            density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
        else
            density = hist2d(data{1},resX, resY,Xrange,Yrange);
        end

        %display for ROI definition
        figure;hIm = imshow(density,'XData',Xrange,'YData',Yrange); axis equal tight off;
        %mask = roipoly;


        %if isempty(roicoords);
        h = imrect;
        roicoords = getPosition(h);
    end              
%     else %h = imrect; %(,roicoords); 
%         roicoords
%         setPosition(h,roicoords);
%     end 
      
%     mask = createMask(h,hIm);
%     params = {};
    corrData = repmat(struct('radius',[],'L',[]),numChannels,1);
    %L = {};
    for ii = 1:numChannels
        [corrData(ii).L,vq,selectedRadius] = Ripley(roicoords,data{ii},maxrad1*nmpixSize,ripleyradius,frame);
        %fname = sprintf('clustermap0%d.tif',ii);
        %imwrite(uint16(vq),fname,'tif','compression','lzw')
    end

function [minX, maxX, minY, maxY] = getBounds(data, Xcol, Ycol)
    minX = min(data(:,Xcol));
    maxX = max(data(:,Xcol));
    minY = min(data(:,Ycol));
    maxY = max(data(:,Ycol));
      
function [L,vq,distnew] = Ripley(roicoords,data,radius,ripleyradius,frame)
    [locs, range, box] = getCoords(roicoords,data);
    maxrad = radius;
    dr = maxrad / 200;
    r = 0:dr:maxrad;
    [K,L] = RipleysK(locs,r,box,0);
    
    if frame == 1
        figure;plot(r,L)
        [distnew,~] = ginput(1);
%         if isempty(ripleyradius);
%             [distnew,~] = ginput(1);                %crosshairs - remove once radius set
%         else
%             distnew = ripleyradius;
%         end
    else
        distnew = ripleyradius;
    end
    %[K,l] = RipleysKPerPoint(locs,distnew,box,0);
    %Lthresh = max(l);
    N = localisation_density(locs,distnew);   
    [ripX,ripY] = meshgrid(range{1},range{2});
    vq = griddata(locs(:,1),locs(:,2),N,ripX,ripY,'v4');     
    if frame == 1
        figure;scatter(locs(:,1),locs(:,2),20,N);

%         [ripX,ripY] = meshgrid(range{1},range{2});
%         vq = griddata(locs(:,1),locs(:,2),N,ripX,ripY,'v4');       
        figure;imagesc(vq); axis xy
    end
        
    
function [coords, range, box] = getCoords(roicoords,data)
%     xyLim = getPosition(h);
%     xmin = xyLim(1);
%     ymin = xyLim(2);
%     width = xyLim(3);
%     height = xyLim(4);
    xmin = roicoords(1);
    ymin = roicoords(2);
%     xmax = xmin + width;
%     ymax = ymin + height;
    xmax = xmin + roicoords(3);
    ymax = ymin + roicoords(4);
    
    dx = 10; %nm
    dy = 10; %nm
    Xrange = xmin:dx:xmax;
    Yrange = ymin:dy:ymax;
    range = {};
    range{1} = Xrange;
    range{2} = Yrange;
    XPosition = data(:,1);
    YPosition = data(:,2);
    isROI = XPosition>xmin & XPosition<xmax & YPosition>ymin & YPosition<ymax;
    X = XPosition(isROI);
    Y = YPosition(isROI);
    coords = [X Y];
    box = [xmin, xmax, ymin, ymax];

    
function N = localisation_density(locs, r)
    DIST = createDistanceMatrix(locs,locs);
    DIST = sort(DIST);
    N = zeros(size(locs,1),1);
    for j = 1: size(locs,1)
        N(j) = (length(find(DIST(2:end,j)<r)));
    end
