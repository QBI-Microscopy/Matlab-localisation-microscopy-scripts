function [Kpp,Lpp] = ripleykperpoint(dataXY,xK,box,method)
% KFUNCTION calculates Ripleys K function
% K = kfunction(dataXY,xK,box,method) - returns vector K containing value
% of Ripley's K-function of dataXY in the distances in xK.
% dataXY - N-by-2 vector where N is number of datapoints. Each row
% corresponds to x and y coordinates of each datapoint
% xK - corresponds to the distances where K function should be computed.
% K is the same size as xK...
% box - rectangular boudnary of the data: box = [xlim1, xlim2, ylim1,
% ylim2]
% method - switch between edge correction. If method=0, no edge correction
% is applied. If method=1, datapoint is used for estimation of K(h) only if
% it is at least h units away from the box

if nargin<4 method=1; end
[N,k] = size(dataXY);
if k~=2 error('dataXY must have two columns'); end

if size(box,1) > 1
    rbox = zeros(N,1);
    for i = 1: N            
        rbox(i) = abs(p_poly_dist(dataXY(i,1),dataXY(i,2), box(:,1),box(:,2)));
    end
    A = polyarea(box(:,1),box(:,2));
else
    rbox = min([ dataXY(:,1)'-box(1);
                    box(2)-dataXY(:,1)';
                    dataXY(:,2)'-box(3);
                    box(4)-dataXY(:,2)']);
    A = (box(2)-box(1))*(box(4)-box(3));
end

%DIST = createDistanceMatrix(dataXY,dataXY);
DIST = pdist2(dataXY,dataXY,'euclidean');
DIST = sort(DIST);

if method == 0 % no correction...
    Kpp = zeros(N,1);
    for k = 1:length(Kpp)
        Kpp(k) = A*sum(DIST(2:end,k)<xK)/N;
        Lpp(k) = sqrt(Kpp(k)/pi) - xK;
    end    
elseif method == 1 % edge correction
    Kpp = zeros(N,1);
    for k=1:length(Kpp)
        I = find(rbox(k)>xK, 1);
        if ~isempty(I)
            Kpp(k) = A*(sum(DIST(2:end,k)<xK))/(length(I));
            Lpp(k) = sqrt(Kpp(k)/pi) - xK;
        end
    end
elseif method == 2 % global edge correction
    W = box(2)-box(1);
    H = box(4)-box(3);
    edge  = 1 - (4/3/pi)*(xK/H + xK/W) + (11/3/pi - 1)*(xK^2/H/W);
    Kpp = zeros(N,1);
    for k=1:length(Kpp)
        Kpp(k) = A*(sum(DIST(2:end,k)<xK))/N/edge;
        Lpp(k) = sqrt(Kpp(k)/pi);
    end  
    
elseif method == 3
    weight = ones(size(DIST));
    for ii = 2: N
        for jj = 1: N
            if DIST(ii,jj) > rbox(ii)
                weight(ii,jj) = 1 - (acos(rbox(ii) / DIST(ii,jj))) / pi;
            end
        end
    end
    Kpp = zeros(N,1);
    Lpp = zeros(N,1);
    for k=1:length(Kpp)
        Kpp(k) = A*(sum((DIST(2:end,k)<xK)./weight(2:end,k)))/N;
    end 
end

