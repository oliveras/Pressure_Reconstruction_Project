    function output = solveFunc(X,x,g,alpha,sigma,h)
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% output = solveFunc(X,x,g,dx,alpha,sigma,h)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Function containing the appropriate equations for solving the nonlocal
% problem.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 2/5/2009


eta = X(1:end-1); etaX = d(eta); etaXX = d(etaX);
M = length(eta); dx = 1/M*2*pi;
N = (M-1)/2;
NVect = (-N:N)';
c = X(end);

[cosEq sinEq] = getSpectralEquations(eta, c, x,g,alpha,sigma,h);

if sigma ==0
    surfaceTensionTerm = 0;
else
    surfaceTensionTerm =  sigma * etaXX./(sqrt( 1+etaX.^2).^3);
end

finalEqn =  conservedQuantity(eta,c,g,surfaceTensionTerm,dx)-alpha;

output = [cosEq;sinEq; finalEqn];