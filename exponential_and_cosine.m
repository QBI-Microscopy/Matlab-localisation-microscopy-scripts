function y = exponential_and_cosine(P, r)
xi = P(1);
A = P(2);
r0 = P(3);
if length(P)<4, C = 1; else C = P(4); end

y = C + A*exp(-r/xi).*cos((pi*r)/2/r0);