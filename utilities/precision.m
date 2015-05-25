function sigma = precision(s,a,N,b)
%   s is the standard deviation of the PSF in nm. This was given in pixels in
%   input and now needs to be converted to nm
%   a is pixel size in image
%   N is the total number of photons measured from the molecule
%   b is the number of photons collected within the fitting window
    snm = s*a;
    sigma = sqrt((snm^2 + a^2/12)/N + (4*sqrt(pi)*snm^3*b^2)/(a*N^2));
end