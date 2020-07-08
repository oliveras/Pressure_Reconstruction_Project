function [eta c finalError] = findSolution(iniGuessEta, iniGuessC, g, h, L, sigma, options)
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [eta c finalError] = findSolution(iniGuess, period, c,sigma,h)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Using an initial guess, the function calls sets parameters for the
% nonlinear solver.  The parameters 'sigma' and 'h' are optional parameters
% and are set to 0 and infinity by default (the case of no surface tension
% and infinite depth.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 2/5/2009


M = length(iniGuessEta);
N = (M-1)/2;
x = (0:(M-1))'/M*L;


eta = iniGuessEta;
c = iniGuessC;
etaX = d(eta);
etaXX = d(etaX);

if sigma ==0
    surfaceTensionTerm = 0;
else
    surfaceTensionTerm =  sigma * etaXX./(sqrt( 1+etaX.^2).^3);
end

% Initial Data and Guess

dx = x(2);


alpha = conservedQuantity(iniGuessEta,iniGuessC,g,surfaceTensionTerm,dx);


initialError  = checkError(iniGuessEta,iniGuessC,sigma,h,g,alpha,surfaceTensionTerm);

if initialError>1e-15
    
    [X] = fsolve(@(X) solveFunc(X,x,g,alpha,sigma,h), [iniGuessEta; iniGuessC],options);
    eta = X(1:end-1);
    c = X(end);
    finalError = checkError(eta,c,sigma,h,g,alpha,surfaceTensionTerm);
else
    X = [eta;c];
    finalError = initialError;
end


