function [x etaNoise etaFilter etaHydro etaKatieH etaFull] = getApproxsFromData(pFilter, c, h, g, L,etaWindow,etaFilter)

M = length(pFilter); N = M/2;
kv = (-N:N-1)'*2*pi/L;
options = optimset('Display','iter','TolFun',1e-15,'TolX',1e-15,'Display','off');

pHat = fftshift(fft(pFilter));


etaHydro = pFilter;


etaKatieH = ifft(ifftshift(cosh(kv*h).*pHat))./(1 - ifft(ifftshift(kv.*pHat.*sinh(kv*h))));


etaFull = real(fsolve(@(etaT) reconstructSurfaceData(etaT,pFilter,c,h,g,L), etaKatieH, options));



interpSize = 10001;
x = linspace(0,L,interpSize);

etaNoise = interpft(etaWindow,interpSize);
etaFilter = interpft(etaFilter,interpSize);
etaHydro = zeroMode(interpft(etaHydro,interpSize)); 
etaKatieH = interpft(etaKatieH,interpSize);
etaFull = interpft(etaFull,interpSize);

