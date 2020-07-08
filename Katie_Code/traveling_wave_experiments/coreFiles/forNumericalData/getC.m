function [c cV] = getC(eta,h,g)

M = length(eta); N = (M-1)/2; x = (0:2*N)'/M*2*pi;


fHat = I(eta,h);
den = 0;
for k = -N:N
    kInd = N+k+1;    
    den = den + 2*exp(i*k*x).*fHat(kInd).*cosh(k*(eta+h));
    for l = -N:N
        lInd = N+l+1;
        den = den + exp(i*(k+l)*x).*fHat(kInd).*fHat(lInd).*exp((k-l)*(eta + h));
    end
end

c2 = -2*g*eta./den;
c = hat(real(sqrt(c2)),0);
cV = real(sqrt(c2));