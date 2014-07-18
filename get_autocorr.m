function [varargout] = get_autocorr(I1 , mask, rmax, flag)
% function [G, r, g, dg, mask] = get_autocorr(I1 , mask, rmax, flag)
% calculates autocorrelation function for two dimensional images
%
% INPUTS
% I1 = image to be autocorrelated
% mask = region of interest, If none or [] specificed, user will be asked
%   to define.  To autocorreate the entire images, use mask = ones(size(I1))
% rmax = maximum r value to correlate in units of pixels. default is 100;
% flag = display flag.  insert 1 to display errorbar(r, g, dg) after
%   computation.
%
% OUTPUTS
% G = two dimensional correlation function.  x and y values range between
%    -rmax:rmax
% r = radius values
% g = angularly averaged autocorrelation function.
% dg = errors on angularly averaged g
% mask = masked used for calculation
%
% NOTE: G(r=0) is just the dot product of the image.  For display purposes,
% G(r=0) is set to zero in the 2D autocorrelation output.  g(r=0) [g(1)]
% retains the proper value.
%
% Last updated 01.26.10 by Sarah Veatch.



if nargin<4, flag = 0; end  % flag for display
if (nargin<3 || isempty(rmax)), rmax=100; end  % distance of maximum correlation recorded
if (nargin<2 || isempty(mask)),    %% draw a mask if needed
    figure;hIm = imshow(I1); axis equal tight off;
    %mask = roipoly;
    h = imrect;
    mask = createMask(h,hIm);
end

%mask = ones(size(I1));
size(mask)
size(I1)
N = sum(sum(I1.*mask));  % number of particles within mask
A = sum(sum(mask));      % area of mask
I1 = double(I1);         % convert to double
L1 = size(I1, 1)+rmax; % size of fft2 (for zero padding)
L2 = size(I1, 2)+rmax; % size of fft2 (for zero padding)

NP = real(fftshift(ifft2(abs(fft2(mask, L1, L2)).^2))); % Normalization for correct boundary conditions
figure;imshow(NP)
G1 = A^2/N^2*real(fftshift(ifft2(abs(fft2(I1.*mask,L1, L2)).^2)))./NP; % 2D G(r) with proper normalization
G = imcrop(G1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);  %only return valid part of G


xvals = ones(1, 2*rmax+1)'*(-rmax:rmax);    %map to x positions with center x=0
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1);    %map to y positions with center y=0
zvals = G;

[theta,r,v] = cart2pol(xvals,yvals, zvals);  % convert x, y to polar coordinates

Ar = reshape(r,1, (2*rmax+1)^2);
Avals = reshape(v,1, (2*rmax+1)^2);
[rr,ind] = sort(Ar);                         % sort by r values
vv = Avals(ind);                             % reindex g
r = 0:floor(max(rr));
max(rr)% the radii you want to extract
[n bin] = histc(rr, r-.5);                   % bin by radius
for j = 1:rmax+1;                            % now get averages
    
    m = bin==j;
    n2 = sum(m);                             % the number of pixels in that bin
    if n2==0, vals(j)=0; er(j)=0;            % if no bins, no data
    else
        g(j) = sum(m.*vv)/n2;               % the average G values in this bin
        dg(j) = sqrt(sum(m.*(vv-g(j)).^2))/n2; % the variance of the mean
    end
end

r = 0:rmax;

%end

G(rmax+1, rmax+1) = 0;

if flag,
    r = 0:rmax;
    figure;errorbar(r(2:length(r)), g(2:length(r)), dg(2:length(r)));
    axis tight
end
varargout{1} = G;
varargout{2} = r;
varargout{3} = g;
varargout{4} = dg;
varargout{5} = mask;
assignin('base','G',G);
assignin('base','r',r);
assignin('base','g',g);
assignin('base','dg',dg);
assignin('base','mask',mask);
