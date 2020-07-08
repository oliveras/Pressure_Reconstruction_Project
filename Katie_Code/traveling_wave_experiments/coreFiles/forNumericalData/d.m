function du = d(u,n);
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
% du/dx = d(u)
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Returns the derivative with respect to x
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 2/5/2009

M = length(u);
N = (M-1)/2;
if nargin<2
    n = 1;
end

%Filter = getFilter(N);
NVect = (-N:N)';
uHat = hat(u);

uXHat = (i*NVect).^n.*uHat;

du = real(invHat(uXHat));
