clear all
close all
clc
format compact

%e01-e05     h=5.05 cm.
%e06-e07     h=3.55 cm.
%e08-e10     h=4.10 cm.

% - - - - - - - - - - - - - - - - - - - - - - - - - 
data = importdata('../../data/e1-9may11-black.txt');
h = 5.05;
g = 980;

p = data(:,1); pMaxInd = find(p==max(p));
eta = data(:,2); etaMaxInd = find(eta==max(eta));

dt = .002;


subplot(2,1,1)
plot(p); hold on
ylabel('Pressure')
xlim(pMaxInd+[-200 200])

subplot(2,1,2)
plot(eta); hold on
ylabel('Surface')

xlim(etaMaxInd+[-200 200])

pWindowInd = pMaxInd+[-200 200];
etaWindowInd = etaMaxInd+[-200 200];



pWindow = p(pWindowInd(1):pWindowInd(2));
etaWindow = eta(etaWindowInd(1):etaWindowInd(2));

pFilter = interpft(pWindow,8);


interpSize = 128;
subplot(2,1,1)
plot(linspace(pWindowInd(1),pWindowInd(2),interpSize),interpft(pFilter,interpSize),'r','LineWidth',2)

