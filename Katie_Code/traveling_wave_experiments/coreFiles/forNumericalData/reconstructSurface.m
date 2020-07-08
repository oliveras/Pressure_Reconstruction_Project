function ans = reconstructSurface(eta,pHat,c,h,g)

M = length(eta); N = (M-1)/2; x = (0:2*N)'/M*2*pi;
kv = (-N:N)'; etaX = d(eta);



for k = -N:N
    j = k + N + 1;
    coshTerm = cosh(kv*(eta(j) + h));
    rhs(j,1) = sum(exp(i*kv*x(j)).*pHat.*coshTerm);
end

lhs = c - sqrt(c^2 - 2*g*eta)./sqrt(1+etaX.^2);
rhs = rhs;
ans = real(lhs-rhs);