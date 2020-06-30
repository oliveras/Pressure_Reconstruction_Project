clear all
close all
clc

% --------------------
% Ploting options string
axes_options = {'interpreter','latex','fontsize',16};
title_options = {'interpreter','latex','fontsize',20};

% --------------------
% Variables in Physical Space
Nx = 100;                       % number of x points
x = linspace(0,2*pi,Nx);        % x grid points
dx = x(2) - x(1);               % determine grid spacing
L = x(end) - x(1);              % determine the period

y = 4*cos(2*x).^3;              % define the periodic function

% --------------------
% Fourier Parameters
Nn = 10;                        % number of Fourier Modes
N = (-Nn:Nn)';                  % a vector of modes

% --------------------
% Calculate the Fourier Coefficients 
intTerm = exp(-1i*N*x).*y;      % create a matrix mesh of the Fourier integrand    
yHat = sum(intTerm(:,2:end),2)*dx/L';   % This is the trapz. rule for periodic functions



% --------------------
% Ignore the imaginary part if the norm is less than 10^(-12)
if norm(imag(yHat))<1e-12
    yHat = real(yHat);
else
    disp('Complex Valued Fourier Coefficients, modify plots accordingly')
end

% --------------------
% Calculate the Fourier Series Representation

sumTerm = exp(1i*N*x).*yHat;
ySum = sum(sumTerm,1);

% --------------------
% Ignore the imaginary part if the norm is less than 10^(-12)
if norm(imag(ySum))<1e-12
    ySum = real(ySum);
else
    disp('Complex Valued Reconstructed y(x), modify plots accordingly')
end




% --------------------
% Plot both the function and the Fourier coefficients
subplot(2,1,1)
plot(x,y,'.-'); hold on
plot(x,ySum)
legend('Original','Reconstructed')
xlabel('$x$',axes_options{:})
ylabel('$y$',axes_options{:})
title('Original function $y(x) = 4\cos^3(2x)$',title_options{:})
xlim([0,L])
grid on


subplot(2,1,2)
stem(N,yHat)
xlabel('$N$',axes_options{:})
ylabel('$\hat{y}$',axes_options{:})
title('Real Part of Fourier Coefficients $\hat{y}$',title_options{:})
xlim([-Nn,Nn])