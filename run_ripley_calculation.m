function [corrData,ZI,distnew,roi_extent] = run_ripley_calculation(varargin)
    
    % extract inputs
    channel1 = varargin{1};
    channel2 = varargin{2};
    resX = varargin{3};
    resY = varargin{4};
    nmPixSize = varargin{5};
    range = varargin{6};
    Xrange = [0 range(1)];
    Yrange = [0 range(2)];
    maxrad = varargin{7};
    ripleyradius = varargin{8};
    roi_extent = varargin{9};
    frame = varargin{10};
    
    % organise input data into cell array
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
        % if there are 2 channels do we want to channel align?
        if ~isempty(t_params{1}) && t_params{2} == 0
            tformData = load(t_params{1});
            TFORM = tformData.TFORM;
            data{1} = tformfwd(TFORM,channel1);
        elseif isempty(t_params{1}) && t_params{2} == 0
            data{1} = channel1;
        elseif isempty(t_params{1}) && t_params{2} == 1
            [in_points,base_points,~,TFORM] = transformChannels(nmPixSize,density);
            data{1} = tformfwd(TFORM,channel1);
        end
        data{2} = channel2;
    end
    
    % if we are on frame 0 make a 2D molecular density map for the 
    % whole dataset - this is used to create an ROI
    if frame == 0
        if numChannels == 2
            density(:,:,1) = hist2d(data{1},resX, resY,Xrange,Yrange);
            density(:,:,2) = hist2d(data{2},resX, resY,Xrange,Yrange);
            density(:,:,3) = zeros(size(density(:,:,1),1),size(density(:,:,1),2));
            dens_min = min(min(min(density(:,:,1))),min(min(density(:,:,2))));
            dens_max = max(max(max(density(:,:,1))),max(max(density(:,:,2))));
        else
            density = hist2d(data{1},resX, resY,Xrange,Yrange);
            dens_min = min(min(density));
            dens_max = max(max(density));
        end

        figure;hIm = imshow(density,'XData',Xrange,'YData',Yrange,'DisplayRange',[0 1]); axis equal tight off;
        h = imrect;
        roi_extent = getPosition(h);
    end
    
    % we are going to make molecular density maps for the ROI
    % and write these as tiff
    roi_xrange = [roi_extent(1) roi_extent(1)+roi_extent(3)];
    roi_yrange = [roi_extent(2) roi_extent(2)+roi_extent(4)];
    roi_xres = ceil((roi_xrange(2) - roi_xrange(1))/nmPixSize);
    roi_yres = ceil((roi_yrange(2) - roi_yrange(1))/nmPixSize);
    [roicoords_ch1, ~] = get_coords_in_roi(roi_extent,data{1});
    if numChannels == 2
        [roicoords_ch2, ~, ~] = get_coords_in_roi(roi_extent,data{2});
        roi_density(:,:,1) = hist2d(roicoords_ch1,roi_xres, roi_yres,roi_xrange,roi_yrange);
        roi_density(:,:,2) = hist2d(roicoords_ch2,roi_xres, roi_yres,roi_xrange,roi_yrange);
        roi_density(:,:,3) = zeros(size(roi_density(:,:,1),1),size(roi_density(:,:,1),2));
    else
        roi_density = hist2d(roicoords_ch1,roi_xres, roi_yres,roi_xrange,roi_yrange);
    end
    if frame > 0
        fname = sprintf('roi_2dhistogram0%d.tif',frame);
        imwrite(uint8(roi_density),fname,'tif','compression','lzw')
        corrData = repmat(struct('radius',[],'L',[]),numChannels,1);
    end
    
    % now do the ripley function calculation
    for ii = 1:numChannels
        % get the localisations within the ROI
        [locs, box] = get_coords_in_roi(roi_extent,data{ii});
        dr = maxrad / 200;
        r = 0:dr:maxrad;
        % calculate k and l functions
        [K,L] = ripleykfunction(locs,r,box,0);
        corrData(ii).L = L;
        % if we are on frame 0 we are doing the calculation for the entire
        % dataset - this can be very slow but is used to get the radius
        % at which there is most clustering
        if frame == 0
            figure;plot(r,L)
            [distnew,~] = ginput(1);
        else
            distnew = ripleyradius;
        end
        % now determine the localisation density in the ROI
        N = ripleykperpoint(locs,distnew,box,0);
        pixelsX = ceil((box(2)-box(1))/50);
        pixelsY = ceil((box(4)-box(3))/50);
        tx = linspace(box(1),box(2),pixelsX);
        ty = linspace(box(3),box(4),pixelsY);
        [XI,YI] = meshgrid(tx,ty);
        ZI = biharmonic_spline_interp2(locs(:,1),locs(:,2),N,XI,YI);   
        if frame == 1
            figure;scatter(locs(:,1),locs(:,2),20,N);      
            figure;imagesc(ZI); axis xy
        end
    end

