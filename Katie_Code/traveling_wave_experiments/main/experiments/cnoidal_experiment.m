clear all
close all
clc

figure
format compact
addpath ../../coreFiles -end


str_options = {'interpreter','latex','fontsize',14};
title_options = {'interpreter','latex','fontsize',18};


%e01-e05     h=5.05 cm.
%e06-e07     h=3.55 cm.
%e08-e10     h=4.10 cm.

% - - - - - - - - - - - - - - - - - - - - - - - - -
data = importdata('../../data/cnshifted-tacp-21aug13-p59-3.txt');

x = data(:,1);
p = data(:,3); pMaxInd = round(mean(find(p==max(p))));
eta = data(:,2); etaMaxInd = round(mean(find(eta==max(eta))));





% - - - - - - - - - - - - - - - - - - - - - - - - -
% Initial Set-Up
% - - - - - - - - - - - - - - - - - - - - - - - - -


pWindowInd = pMaxInd+[-250 250];
etaWindowInd = etaMaxInd+[-250 250]-4;

%2292 2607;
pWindowIndVals = [2442 2787; 1971+173 2606+173; 2609+173 3550+173];

pWindowIndVals = [2000 length(eta)];
for scales = 1;
    for windowIter = 1
        pWindowInd = pWindowIndVals(windowIter,:);
        
        
        etaWindowInd = pWindowInd;
        pWindow = zeroMode(p(pWindowInd(1):pWindowInd(2)));
        etaWindow = zeroMode(eta(etaWindowInd(1):etaWindowInd(2)));
        
        
        if windowIter == 1
            pFilter = zeroMode(filterData(pWindow,32,200));
            etaFilter = zeroMode(filterData(etaWindow,64,200));
        elseif windowIter ==2
            pFilter = zeroMode(filterData(pWindow,10,18));
            etaFilter = zeroMode(filterData(etaWindow,10,18));
        else
            pFilter = zeroMode(filterData(pWindow,20,64));
            etaFilter = zeroMode(filterData(etaWindow,20,30));
        end
        
        
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - -
        
        % - - - - - - - - - - - - - - - - - - - - - - - - -
        % Reconstruct the Surface
        % - - - - - - - - - - - - - - - - - - - - - - - - -
        
        g = 981; h = 6.170584700000010;
        
        
        %%%% Change this to the measured speed.
        
        c0 = sqrt(g*h); dt = .005;
        %c0 = sqrt(g*h)+55.523961849884493;
        etaT = interpft(etaFilter,length(etaWindow));
        a = max(etaT)-min(etaT); epsilon = a/h;
        
        c = c0*(1 + epsilon/2);
        L = length(etaWindow)*dt*c;
        K = 2*pi/L;
        cN = sqrt(g*tanh(K*h)/K);
        c = .95*cN;
        N = length(pFilter)/2;
        c = cN*1.5;
        c = sqrt(g*h)*scales;
        c = c0;
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        [x etaNoise etaFilter etaHydro etaKatieH etaFull] = getApproxsFromData(pFilter, c, h, g, L,etaWindow,etaFilter);
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
        
        
        
    end
end

%legend('True .5','Heuristic','Nonlocal','True 1','Heuristic','Nonlocal','True 2','Heuristic','Nonlocal')

%save(fileName)

