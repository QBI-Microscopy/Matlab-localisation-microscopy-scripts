function y = exponential(P, r)
xi = P(1);
A = P(2);
if length(P)<3, C = 0; else C = P(3); end

y = A*exp(-r/xi)+C;