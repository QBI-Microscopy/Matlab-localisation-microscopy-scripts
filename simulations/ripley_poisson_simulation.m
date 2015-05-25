rmax = 1000;
N = 200;
dr = rmax / (N-1);
radius = 0:dr:rmax;
Nobjects = 1000;
numSimulations = 1;    
Lrand = zeros(numSimulations,length(radius));
for r = 1: numSimulations
    rand_centers = [];
    randX = random('unif',box(1),box(2),[Nobjects,1]);
    randY = random('unif',box(3),box(4),[Nobjects,1]);

    rand_centers = [randX,randY];
    if r == 1
        figure;plot(rand_centers(:,1),rand_centers(:,2),'r+')
    end

    xmin = min(rand_centers(:,1));
    ymin = min(rand_centers(:,2));
    xmax = max(rand_centers(:,1));
    ymax = max(rand_centers(:,2));
    box = [xmin, xmax, ymin, ymax];
    [K,Lrand(r,:)] = ripleykfunction(rand_centers,radius,box,1);
end
figure
plot(radius,Lrand(1,:),'k-')
