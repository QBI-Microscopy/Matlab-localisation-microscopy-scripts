clc, clear all

%% parameters
sizeX = 500; %region size in arbitrary units
sizeY = 500;
density = 1e-3;
N = density*sizeX*sizeY; %particle number
NumRClusters = 25;
clusterFact = 1; %change the relative number of red and green clusters (vary between 0 -> 1)
NumGClusters = ceil(clusterFact*NumRClusters);
NParticlesPerCluster = floor(N/NumRClusters);
rRadius = 10;
radiusFact = 1; %change the relative size of red to green clusters (vary between 0 -> 1)
gRadius = ceil(radiusFact*rRadius);

%% Generate centers
centers = [];
for ii = 1: NumRClusters
    new_centers = [random('unif',0,sizeX),random('unif',0,sizeY)];
    centers = [centers;new_centers];
end
%% Generate clusters
red_centers = [];
for ii = 1: NumRClusters
    [r_x,r_y] = circ_rand(NParticlesPerCluster,rRadius);
    r_uniform_centers = [(r_x + centers(ii,1)),(r_y + centers(ii,2))];
    red_centers = [red_centers;r_uniform_centers];
end
green_centers = [];
for ii = 1: NumGClusters
    [g_x,g_y] = circ_rand(NParticlesPerCluster,gRadius);
    g_uniform_centers = [(g_x + centers(ii,1)),(g_y + centers(ii,2))];
    green_centers = [green_centers;g_uniform_centers];
end
figure
hold on
plot(red_centers(:,1),red_centers(:,2),'r+');
plot(green_centers(:,1),green_centers(:,2),'g+');
hold off


%% Create images
r_density_centers = hist2d(red_centers,sizeX, sizeY,[0 sizeX],[0 sizeY]);
g_density_centers = hist2d(green_centers,sizeX, sizeY,[0 sizeX],[0 sizeY]);

%% calculate the correlation and the fit
maxrad1 = 200;
Imsize1=min(size(r_density_centers)); 
if Imsize1<1.25*maxrad1
    maxrad1=round(Imsize1/1.25);
end
mask = ones(sizeY,sizeX);
[C, r_centers, c, dc, maskout] = get_crosscorr(r_density_centers, g_density_centers,mask, maxrad1, 1);

