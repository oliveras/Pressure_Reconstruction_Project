clear all
close all
clc
format long
format compact
addpath ../../coreFiles -end
load ../output/h_01_N_64_etaTrue_Pressure_06_29_2011


ampInd = 100;

p = pOut(:,ampInd); M = length(p);
cVals = cTrue(ampInd)+erf((0:10)/10);
%cVals = erfc(linspace(-2.5,2.5,10))*;
etaT = etaTrue(:,ampInd);
x = (0:(M-1))'/M*2*pi;
for j = 1:length(cVals)
    c = cVals(j);
    eta = real(fsolve(@(eta) reconstructSurface(eta,hat(c-sqrt(c^2-2*p)),c,h,g), p, options));
    errorVal(j) = norm(eta - etaT);
    plot(x,[eta etaT])
    pause(.3)
end

