clear all
close all
clc


format compact
addpath ../../coreFiles -end

str_options = {'interpreter','latex','fontsize',14};
title_options = {'interpreter','latex','fontsize',18};


%e01-e05     h=5.05 cm.
%e06-e07     h=3.55 cm.
%e08-e10     h=4.10 cm.
% - - - - - - - - - - - - - - - - - - - - - - - - - 
data = importdata('../../data/e1-9may11-black.txt');

p = data(:,1); pMaxInd = round(mean(find(p==max(p))));
eta = data(:,2); etaMaxInd = round(mean(find(eta==max(eta))));

t = 0:.002:30; t = t(1:end-1);


% - - - - - - - - - - - - - - - - - - - - - - - - - 
% Let's plot some data
% - - - - - - - - - - - - - - - - - - - - - - - - - 

% Plot of raw pressure and surface data
subplot(2,1,1)
plot(t,p); hold on; plot(t,eta); hold off;
legend('Pressure Data - Raw','Surface Data - Raw','location','best')
xlabel('Time',str_options{:})
ylabel('Measured Value in cm',str_options{:})
title('Raw Pressure and Surface Data',title_options{:})
grid on


subplot(2,1,2)

M = length(p);
if mod(M,2)==0
    N = M/2;
    Nv = -N:(N-1);
else
    N = (M-1)/2;
    Nv = -N:N;
end
Nv = -7500:7499;
filter = abs(Nv<20);
pHat = fftshift(fft(p)).*filter;
stem(Nv,real(pHat)/M)
hold on
stem(Nv,imag(pHat)/M)
hold off
leg = legend('$\mathcal{R}\lbrace\hat{p}\rbrace$','$\mathcal{I}\lbrace\hat{p}\rbrace$');
set(leg,'location','best',str_options{:})
xlabel('Index',str_options{:})
ylabel('Value',str_options{:})
title('FFT of the Raw Pressure Data',title_options{:})
set(gcf,'position',[446         283        1779         949])

return

% - - - - - - - - - - - - - - - - - - - - - - - - - 
% Initial Set-Up
% - - - - - - - - - - - - - - - - - - - - - - - - - 
pWindowInd = pMaxInd+[-250 250];
etaWindowInd = etaMaxInd+[-250 250]-4;

pWindow = zeroMode(p(pWindowInd(1):pWindowInd(2)));
etaWindow = zeroMode(eta(etaWindowInd(1):etaWindowInd(2)));


pFilter = zeroMode(filterData(pWindow,4,10));
etaFilter = zeroMode(filterData(etaWindow,6,12));

% - - - - - - - - - - - - - - - - - - - - - - - - - 
expNum = 1;
fileName = ['./output/e',sprintf('%02d',expNum),'_output'];


% - - - - - - - - - - - - - - - - - - - - - - - - - 
% Reconstruct the Surface
% - - - - - - - - - - - - - - - - - - - - - - - - - 
hVals = [5.05 5.05 5.05 5.05 5.05 3.55 3.55 4.1 4.1 4.1];

g = 981; h = hVals(expNum); c0 = sqrt(g*h); dt = .002;

etaT = interpft(etaFilter,length(etaWindow));
a = max(etaT)-min(etaT); epsilon = a/h;

c = c0*(1 + epsilon/2);
L = length(etaWindow)*dt*c;
K = 2*pi/L;
cN = sqrt(g*tanh(K*h)/K);
c = .95*cN;
N = length(pFilter)/2;

% - - - - - - - - - - - - - - - - - - - - - - - - - -
[x etaNoise etaFilter etaHydro etaKatieH etaFull] = getApproxsFromData(pFilter, c, h, g, L,etaWindow,etaFilter);

subplot(2,1,2)
colors = parula(10);

plot(x,etaFilter,'k','LineWidth',2); hold on
plot(x,etaHydro)%,'Color',colors(1,:),'MarkerFaceColor',colors(1,:))
plot(x,etaKatieH)%,'Color',colors(4,:),'MarkerFaceColor',colors(4,:))
plot(x,etaFull)%,'Color',colors(85,:),'MarkerFaceColor',colors(8,:))
xlabel('\xi - cm')
ylabel('Surface Elevation \eta - cm')
axis tight

leg = legend('Measured Surface','$p = \rho g h$','Katie - Heuristic','Nonlinear Map','location','best');
set(leg,str_options{:});

set(gcf,'position',[446         283        1779         949])


save(fileName)

