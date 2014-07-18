function [corrData] = run_correlation_and_fit(varargin) %(input1, input2, Xcol, Ycol, res, nm pix size, range, correlation type, fittype, radius)

    channel1 = varargin{1};
    channel2 = varargin{2};
    res = varargin{3}
    nmpixSize = varargin{4};
    t_params = varargin{5};
    range = varargin{6};
    Xrange = [0 range(1)];
    Yrange = [0 range(2)];
    correlation = varargin{7};
    fit = varargin{8};
    maxrad1 = varargin{9};
      
    %calculate histograms
    if ~isempty(channel2)
        density(:,:,1) = hist2d(channel1,res(1), res(2),Xrange,Yrange);
        density(:,:,2) = hist2d(channel2,res(1), res(2),Xrange,Yrange);
        density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
    else
        density = hist2d(channel1,res(1), res(2),Xrange,Yrange);
    end    
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
        %calculate histograms
    if numChannels == 2
        density(:,:,1) = hist2d(data{1},res(1), res(2),Xrange,Yrange);
        density(:,:,2) = hist2d(data{2},res(1), res(2),Xrange,Yrange);
        density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
    else
        density = hist2d(data{1},res(1), res(2),Xrange,Yrange);
    end
    
    %display for ROI definition
    figure;hIm = imshow(density,'XData',Xrange,'YData',Yrange); axis equal tight off;
    %mask = roipoly;
    h = imrect;
    pos = getPosition(h)
    mask = createMask(h,hIm);
    figure;imshow(mask)
    params = {};
    corrData = repmat(struct('twoDcorr',[],'radius',[],'correlation',[],'error',[],'mask',[],'type',[],'L',[]),numChannels,1);
    L = {};
    for ii = 1:numChannels
        [corrData(ii).L,vq] = Ripley(h,data{ii},maxrad1*nmpixSize);
        fname = sprintf('clustermap0%d.tif',ii);
        imwrite(uint16(vq),fname,'tif','compression','lzw')
    end
    % the input parameters for the fit depend on the fit type
    switch fit
        case 'exponential_and_gaussian'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
            P(3) = 20; %sqrt(2)* the PSF of the image. (the half with of the gaussian) in nm
            P(4) = 1; % the surface density of the probe (in 1/um^2)
        case 'exponential_and_cosine'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
            P(3) = 20; %sqrt(2)* the PSF of the image. (the half with of the gaussian) in nm
        case 'exponential'
            P = [];
            P(1) = 800; % the decay of the exponential (in nm)
            P(2) = 10; %the amplitude of the exponential (y intercept)
    end

    %calculate correlation
    switch correlation
        case 'auto'
            for iChan = 1: numChannels
                Imsize1=min(size(density(:,:,1))); 
                if Imsize1<1.25*maxrad1
                    maxrad1=round(Imsize1/1.25);
                end
                [G, r, g, dg, maskout] = get_autocorr(density(:,:,iChan), mask, maxrad1, 1);
                if isnan(g)
                    errordlg('auto-correlation calculation failed. try changing radius','modal');
                    return
                else
                    params{iChan} = fitData(r,g,dg,nmpixSize,fit,P);
                    corrData(iChan).twoDcorr = G;
                    corrData(iChan).radius = r;
                    corrData(iChan).correlation = g;
                    corrData(iChan).error = dg;
                    corrData(iChan).mask = maskout;
                    corrData(iChan).type = 'auto';
                end
            end
        case 'cross'
            Imsize1=min(size(density(:,:,1))); 
            if Imsize1<1.25*maxrad1
                maxrad1=round(Imsize1/1.25);
            end
            
            if numChannels == 2
                [C, r, c, dc, maskout] = get_crosscorr(density(:,:,1), density(:,:,2), mask, maxrad1, 1);
            else
                errordlg('two channels are required to calculate cross-correlation','modal');
                return
            end
            
            if isnan(c)
                errordlg('cross-correlation calculation failed. try changing radius','modal');
                return
            else
                params{1} = fitData(r,c,dc,nmpixSize,fit,P);
                corrData(1).twoDcorr = C;
                corrData(1).radius = r;
                corrData(1).correlation = c;
                corrData(1).error = dc;
                corrData(1).mask = maskout;
                corrData(1).type = 'cross';
            end
    end
   


%fit the data
function params = fitData(x,y,err,pixSize,type,Pin)

    switch type
        case 'exponential_and_gaussian'
            params = repmat(struct('cluster_size',[],'magnitude',[],'density',[],'sigma',[]),1,1);
        case 'exponential_and_cosine'
            params = repmat(struct('cluster_size',[],'magnitude',[],'r0',[]),1,1);
        case 'exponential'
            params = repmat(struct('cluster_size',[],'magnitude',[]),1,1);    
    end

    x = x .* pixSize;
    figure;errorbar(x(2:end), y(2:end), err(2:end), '.')
    switch type
        case 'exponential_and_gaussian'
            P0 = [Pin(1), Pin(2), Pin(3), Pin(4)];

            P = lsqcurvefit('exponential_and_gaussian', P0, x(2:end), y(2:end)-1);

            hold on
            plot(x, exponential_and_gaussian(P, x)+1, 'r')
            hold off

            legend ('data', 'fit')

            params.cluster_size = P(1) %characteristic size of structure
            params.magnitude = P(2) % magnitude of clustering
            params.sigma = P(3)/sqrt(2) %in nm
            params.density = P(4) %in 1/um^
        case 'exponential_and_cosine'
            P0 = [Pin(1), Pin(2), Pin(3)];
            P = lsqcurvefit('exponential_and_cosine', P0, x(2:end), y(2:end)-1);

            hold on
            plot(x, exponential_and_cosine(P, x)+1, 'r')
            hold off

            legend ('data', 'fit')

            params.cluster_size = P(1) %characteristic size of structure
            params.magnitude = P(2) % magnitude of clustering
            params.r0 = P(3) %in nm  
        case 'exponential'
            P0 = [Pin(1), Pin(2)];
            P = lsqcurvefit('exponential', P0, x(2:end), y(2:end)-1);

            hold on
            plot(x, exponential(P, x)+1, 'r')
            hold off

            legend ('data', 'fit')

            params.cluster_size = P(1) %characteristic size of structure
            params.magnitude = P(2) % magnitude of clustering                 
    end


% function [density, n, m, X, Y] = calcHistogram(N, data, dataCols, res)
%     %N = num rows in data set
%     Xcol = dataCols(1);
%     Ycol = dataCols(2);
%     subset = true(N,1);
%     [minX, maxX, minY, maxY] = getBounds(data, Xcol, Ycol);
%     NP = sum(subset);
%     XPosition = getX(data, Xcol);
%     YPosition = getY(data, Ycol);
%     n=linspace(minX,maxX,res); m=linspace(minY,maxY,res);
% 
%     if NP>0
%         RR = round((res-1)*[...
%             (XPosition(subset)-minX)/(maxX-minX) ...
%             (YPosition(subset)-minY)/(maxY-minY) ])+1;
%         RRok = all(RR<=res,2) & all(RR>=1,2) ;
%         pxx=(maxX-minX)/res; pxy=(maxY-minY)/res;
%         pxArea=pxx * pxy;
%         density = accumarray(RR(RRok,:),1,[res,res])/pxArea;
%         X = repmat(n,res,1); Y = repmat(m',1,res);
%     else
%         density = zeros(res+1);
%         X = repmat(n,res+1,1); Y = repmat(m',1,res+1);
%     end
    
function [minX, maxX, minY, maxY] = getBounds(data, Xcol, Ycol)
    minX = min(data(:,Xcol));
    maxX = max(data(:,Xcol));
    minY = min(data(:,Ycol));
    maxY = max(data(:,Ycol));
      
function [L,vq] = Ripley(h,data,radius)
    [locs, range, box] = getCoords(h,data);
    maxrad = radius;
    dr = maxrad / 200;
    r = dr:dr:maxrad;
    [K,L] = RipleysK(locs,r,box,0);
    [Knew,Lnew] = ripleykfunction(locs,r,box,0);
    Lrand = ripleysimulation(100,r,locs,box);
    figure;
    hold on
    plot(r,L,'b-')
    plot(radius,mean(Lrand,1)+2*std(Lrand,1),'r-')
    plot(radius,mean(Lrand,1)-2*std(Lrand,1),'r-')
    hold off
    [distnew,~] = ginput(1);
    N = ripleykperpoint(locs, distnew, box, 0);
    figure;scatter(locs(:,1),locs(:,2),20,N);
    
    [ripX,ripY] = meshgrid(range{1},range{2});
    vq = griddata(locs(:,1),locs(:,2),N,ripX,ripY,'v4');       
    figure;imagesc(vq); axis xy
    
function Lrand = ripleysimulation(numSimulations,radius,coords,box)

    Nobjects = size(coords,1);

    Lrand = zeros(numSimulations,length(radius));
    for r = 1: numSimulations
        randX = random('unif',min(coords(:,1)),max(coords(:,1)),[Nobjects,1]);
        randY = random('unif',min(coords(:,2)),max(coords(:,2)),[Nobjects,1]);
        rand_centers = [randX,randY];
        if ~isempty(rand_centers)
            [~,Lrand(r,:)] = ripleykfunction(rand_centers,radius,box,0);
        else
            Lrand(r,:) = NaN;
        end
    end
    % This gives a 95% C.I.
    % LrandOut{iROI} = Lrand;
    RadiusOut = radius;
    SDPlus = (mean(Lrand,1)+2*std(Lrand,1))';
    SDMinus = (mean(Lrand,1)-2*std(Lrand,1))';
    
function [coords, range, box] = getCoords(h,data)
    xyLim = getPosition(h);
    xmin = xyLim(1);
    ymin = xyLim(2);
    width = xyLim(3);
    height = xyLim(4);
    xmax = xmin + width;
    ymax = ymin + height;
    
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