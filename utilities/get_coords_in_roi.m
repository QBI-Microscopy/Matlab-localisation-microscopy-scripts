function [coords, box] = get_coords_in_roi(roicoords,data)
    xmin = roicoords(1);
    ymin = roicoords(2);
    xmax = xmin + roicoords(3);
    ymax = ymin + roicoords(4);

    XPosition = data(:,1);
    YPosition = data(:,2);
    isROI = XPosition>xmin & XPosition<xmax & YPosition>ymin & YPosition<ymax;
    X = XPosition(isROI);
    Y = YPosition(isROI);
    coords = [X Y];
    box = [xmin, xmax, ymin, ymax];
