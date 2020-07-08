clear all
close all
clc


format compact
addpath ../../coreFiles -end

% - - - - - - - - - - - - - - - - - - - - - - - - - 
buoyNumber = 32411;
%run(strcat('../buoyData/buoy',num2str(buoyNumber),'_data.m'));

load ../buoyData/buoy32411_Dart_Data
pMaxInd = round(mean(find(p==max(p))));



% - - - - - - - - - - - - - - - - - - - - - - - - - 
% Initial Set-Up
% - - - - - - - - - - - - - - - - - - - - - - - - - 

g = 9.81; % meters / s^2
h = 3243; % meters
pressureConversion = 670; % mm / psi

pWindowInd = pMaxInd+[-100 300];
tFilter = t(pWindowInd(1):pWindowInd(2));
%plot(tFilter,p(pWindowInd(1):pWindowInd(2)));


pWindow = zeroMode(p(pWindowInd(1):pWindowInd(2)));
pFilter = zeroMode(interpft(pWindow,64));


plot(tFilter,pWindow,'k','LineWidth',2)
hold on;
plot(tFilter,interpft(pFilter,length(tFilter)),'k-.','LineWidth',2)


c = sqrt(g*h); 
dt = 15;


L = length(tFilter)*dt*c;
K = 2*pi/L;
M = length(pFilter);
N = M/2;

kv = (-N:N-1)'*2*pi/L;



pHat = fftshift(fft(pFilter));
pHatN = fftshift(fft(c-sqrt(c^2-2*pFilter)));


pX = real(ifft(ifftshift(1i*kv.*pHat)));
pXX = real(ifft(ifftshift(-kv.^2.*pHat)));
pXXX = real(ifft(ifftshift(-1i*kv.^3.*pHat)));
pXXXX = real(ifft(ifftshift(kv.^4.*pHat)));
etaLinearHat = cosh(kv*h).*pHat;
etaLinear = ifft(ifftshift(etaLinearHat));

temp = ifft(ifftshift(kv.*pHat.*sinh(kv*h)));

numTemp = ifft(ifftshift(cosh(kv*h).*pHatN));
denTemp = ifft(ifftshift(sinh(kv*h).*kv.*pHatN));
etaKatieH = etaLinear./(1 - temp);

hold on;

plot(tFilter,interpft(etaLinear,length(tFilter)),'b-','LineWidth',1.5)
plot(tFilter,interpft(etaKatieH,length(tFilter)),'r-.','LineWidth',1.5)
axis tight

ylim(1.1*[min(etaKatieH) max(etaKatieH)])
legend('Measured Data','Filtered Data','Transfer Function','Heuristic Formula')

xlabel('Date - Julian')
ylabel('Height - meters')
title('Buoy #32411 - Chile 2010 Tsunami - depth = 3243 meters')
set(gcf,'Position', [324   332   923   375])