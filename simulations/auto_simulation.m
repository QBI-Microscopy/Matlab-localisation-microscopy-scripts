clc, clear all
%% parameters
sizeX = 500; %region size in arbitrary units
sizeY = 500;
density = 1e-3;
N = density*sizeX*sizeY; %num of particles
NumClusters = 25; %number of domains
NParticlesPerCluster = floor(N/NumClusters);
clusterRadius = 25; %radius of domain
psfSigma = 2; %PSF radius

%% define some random point to simulate background
rand_centers = [];
for ii = 1: 100
    new_rand_centers = [random('unif',0,sizeX),random('unif',0,sizeY)];
    rand_centers = [rand_centers;new_rand_centers];
end

%% define domain centers from uniformly distributed random numbers
centers = [];
for ii = 1: NumClusters
    new_centers = [random('unif',0,sizeX),random('unif',0,sizeY)];
    centers = [centers;new_centers];
end

%% define clusters - circular domains of uniformly distributed random numbers
red_centers = [];
for ii = 1: NumClusters
    [r_x,r_y] = circ_rand(NParticlesPerCluster,clusterRadius);
    r_uniform_centers = [(r_x + centers(ii,1)),(r_y + centers(ii,2))];
    red_centers = [red_centers;r_uniform_centers];
end

%% define over-counting - overcounting ratio (OCR) = 1 picked from random numbers with normal distribution
red_signals = red_centers;
for ii = 1: size(red_centers,1)
    gauss_centers1 = [normrnd(red_centers(ii,1),psfSigma),normrnd(red_centers(ii,2),psfSigma)];
    gauss_centers2 = [normrnd(red_centers(ii,1),psfSigma),normrnd(red_centers(ii,2),psfSigma)];
    red_signals = [red_signals;gauss_centers1;gauss_centers2];
end

%% define randomly distributed population - uniform distribtion
green_centers = [sizeX.*rand(N,1),sizeY.*rand(N,1)];

%% over-counting of random distribution - OCR = 1 picked from random numbers with normal distribution
green_signals = green_centers;
for ii = 1: size(red_centers,1)
    gauss_centers1 = [normrnd(green_centers(ii,1),psfSigma),normrnd(green_centers(ii,2),psfSigma)];
    gauss_centers2 = [normrnd(green_centers(ii,1),psfSigma),normrnd(green_centers(ii,2),psfSigma)];
    green_signals = [green_signals;gauss_centers1;gauss_centers2];
end

%% calculate images - 2d histograms
r_density_centers = hist2d(red_centers,sizeX, sizeY,[0 sizeX],[0 sizeY]);
r_density_signals = hist2d(red_signals,sizeX, sizeY,[0 sizeX],[0 sizeY]);
g_density_centers = hist2d(green_centers,sizeX, sizeY,[0 sizeX],[0 sizeY]);
g_density_signals = hist2d(green_signals,sizeX, sizeY,[0 sizeX],[0 sizeY]);
%centers_and_domains = hist2d(rand_and_domains,sizeX, sizeY,[0 sizeX],[0 sizeY]);
%% calculate the correlation and the fit
maxrad1 = 50; %in pixels
Imsize1=min(size(r_density_centers)); 
if Imsize1<1.25*maxrad1
    maxrad1=round(Imsize1/1.25);
end
mask = ones(sizeY,sizeX);
[Grc, r_centers, grc, dg, maskout] = get_autocorr(r_density_centers, mask, maxrad1, 0);

[Grs, r_signals, grs, dg, maskout] = get_autocorr(r_density_signals, mask, maxrad1, 0);

[Ggc, g_centers, ggc, dg, maskout] = get_autocorr(g_density_centers, mask, maxrad1, 0);

[Ggs, g_signals, ggs, dg, maskout] = get_autocorr(g_density_signals, mask, maxrad1, 0);

%[Grands, rand_signals, grands, dg, maskout] = get_autocorr(centers_and_domains, mask, maxrad1, 0);

% replot with the centers and signals overlaid
figure;hold on; 
plot(r_centers(2:end),grc(2:end),'r-'); 
plot(r_signals(2:end),grs(2:end),'b-');
legend('centers','over-counted')
title('Domains');
xlabel('Radius');
ylabel('G(r)');
hold off
figure;hold on; 
plot(g_centers(2:end),ggc(2:end),'r-');
plot(g_signals(2:end),ggs(2:end),'b-');
legend('centers','over-counted')
title('Random distribution');
xlabel('Radius');
ylabel('G(r)');
hold off
% figure;hold on; 
% plot(r_centers(2:end),grc(2:end),'r-');
% plot(rand_signals(2:end),grands(2:end),'b-');
% legend('centers','centers + rand backround')
% title('Centers + Random background');
% xlabel('Radius');
% ylabel('G(r)');
% hold off