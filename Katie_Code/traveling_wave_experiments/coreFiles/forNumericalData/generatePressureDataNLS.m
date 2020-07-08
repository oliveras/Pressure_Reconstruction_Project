function lhs = generatePressureDataNLS(p,eta,etaX,N,h,c)


pHat = hat(sqrt(c^2 - 2*p));

g=1;
kv = (-N:N)';
for k = -N:N
    ind = N + k + 1;
    coshTerm = cosh(kv*eta(ind)) + tanh(kv*eta(ind)).*sinh(kv*h);
    intTerm = invHat(coshTerm.*pHat);
    lhs(ind,1) = sqrt(c^2-2*g*eta(ind))/sqrt(1 + etaX(ind)^2) - intTerm(ind);
end