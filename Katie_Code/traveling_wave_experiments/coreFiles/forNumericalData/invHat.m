function u = invHat(uHat);
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% u = invHat(cosTerm,sinTerm)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% Determines the function 'u' based on the sine and cosine Fourier
% coefficients.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 2/5/2009


i = sqrt(-1);

M = length(uHat);
N = (M-1)/2;kv = (-N:N)';
x = (0:2*N)'/M*2*pi;
u = zeros(size(uHat));

for k=-N:N
    j = k + N + 1;
    u = u + exp(i*kv(j)*x)*uHat(j);
    
end
u = real(u);