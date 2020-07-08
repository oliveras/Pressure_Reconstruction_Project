function uHat = hat(p,n);
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% uHat = hat(u,k)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% Determines the Fourier coefficients of the Fourier series for the 
% periodic function 'u'.  The argument 'k' is optional.  When 'k' is
% defined, uHat only returns the Fourier coefficient corresponding to the
% 'exp(i*k*x)' term.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 2/5/2009

M = length(p); N = (M-1)/2; x = (0:2*N)'/M*2*pi;kv = (-N:N)';

for k = -N:N
    j = k + N + 1;
    pHat(j,1) = integrate(x,exp(-i*k*x).*p);
end

if nargin > 1
    uHat = pHat(n+N+1,1);
else
    uHat = pHat;
end
        
uHat = (uHat)/2/pi;

