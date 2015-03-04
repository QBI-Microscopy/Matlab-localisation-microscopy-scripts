function N = localisation_density(locs, r)
    DIST = createDistanceMatrix(locs,locs);
    DIST = sort(DIST);
    N = zeros(size(locs,1),1);
    for j = 1: size(locs,1)
        N(j) = (length(find(DIST(2:end,j)<r)));
    end

