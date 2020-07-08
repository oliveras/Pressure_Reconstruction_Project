function [cosEq sinEq] = getSpectralEquations(eta,c, x, g, amp, sigma, h)

N = (length(eta)-1)/2;
dx = x(2);
etaX = d(eta);
etaXX = d(etaX);



% - - - - - - - - - - - - - - -
% Surface Tension Term
% - - - - - - - - - - - - - - -

if sigma ==0
    surfaceTensionTerm = 0;
else
    surfaceTensionTerm =  sigma * etaXX./(1 + etaX.^2).^(3/2);
end


% = = = = = = = = = = = = = = =
% Finite Depth
% = = = = = = = = = = = = = = =


for k = 1:N
    sinhTerm = sinh(k*eta) + tanh(k*h)*cosh(k*eta);
    term = sinhTerm.*sqrt((1 + etaX.^2).*(c^2 - 2*g*eta + 2*surfaceTensionTerm)); 
    termHat = integrate(x,exp(i*k*x).*term);
    cosEq(k,1) = real(termHat);
    sinEq(k,1) = imag(termHat);
end


