function ans = reconstructSurfaceData(eta,p,c,h,g,L)

N = length(eta)/2;
pHat = fftshift(fft(c-sqrt(c^2-2*g*p)));


kv = (-N:N-1)'*2*pi/L; 

etaX = real(ifft(ifftshift(1i*kv.*fftshift(fft(eta)))));


for k = -N:N-1
    j = k + N + 1;
    temp = ifft(ifftshift(pHat.*cosh(kv*(eta(j)+h))));
    rhs(j,1) = temp(j);

end


lhs = c-sqrt(c^2 - 2*g*eta)./sqrt(1+etaX.^2);
rhs = rhs;
ans = (lhs-rhs);