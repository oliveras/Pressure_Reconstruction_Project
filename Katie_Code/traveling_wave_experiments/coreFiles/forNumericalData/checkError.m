function error = checkError(eta,c,sigma,h,g,alpha,surfaceTensionTerm)
%
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% error = checkError(eta,c,sigma,h)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Returns the error from the integral equation.  The parameters 'sigma' and
% 'h' are optional.  They are set to '0' and 'infinitiy' to test the case
% with no surface tension and infinite depth.  Gravity 'g' is assumed to be
% 1, and the period 'L' is assumed to be 2*pi.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 2/5/2009



M = length(eta);
N = (M-1)/2;
x = (0:2*N)'/M*2*pi;
dx = x(2);
NVect = (-N:N)';

[cosEq sinEq] = getSpectralEquations(eta, c, x,g,alpha,sigma,h);
finalEqn =  conservedQuantity(eta,c,g,surfaceTensionTerm,dx)-alpha;

error = max(abs([cosEq;sinEq]));