function alpha = conservedQuantity(eta,c,g,surfaceTensionTerm,dx);

%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% alpha = conservedQuantity(eta,c,g,surfaceTensionTerm,dx)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% This creates a conserved quantity to serve as the final constraint in the
% nonlinear solvers.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 2/26/2009
M = length(eta);
N = (M-1)/2; x = (0:2*N)'/M*2*pi;
etaX = d(eta);
etaHat = hat(eta);

alpha(1,1) = max(eta);
alpha(2,1) = integrate(x,eta);

