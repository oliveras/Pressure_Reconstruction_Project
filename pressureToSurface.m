clear all
close all
% clc

format long

% - - - - - - - - - - - - - - - - - - - 
% Model Parameters
% - - - - - - - - - - - - - - - - - - - 
Nx = 100;                       % number of x points
x = linspace(0,2*pi,Nx);        % x grid points
dx = x(2) - x(1);               % determine grid spacing
L = x(end) - x(1);              % determine the period
h = 0.5;                        % depth of the fluid
g = 1;                          % gravity
rho =1;                         % fluid density

% - - - - - - - - - - - - - - - - - - - 
% Fourier Parameters
% - - - - - - - - - - - - - - - - - - - 
Nn = 10;                        % number of Fourier Modes
N = (-Nn:Nn)';                  % a vector of modes

% - - - - - - - - - - - - - - - - - - - 
% Plotting Options
% - - - - - - - - - - - - - - - - - - - 
axes_options = {'interpreter','latex','fontsize',12};

% - - - - - - - - - - - - - - - - - - - 
% Initial Guess for Newton Method
% - - - - - - - - - - - - - - - - - - - 
iniGuessEta = .0001*cos(x);
iniGuessEtaHat = hat(iniGuessEta,x,N);
iniGuessC = sqrt(g*tanh(h));
% iniGuessX = [iniGuessEtaHat;iniGuessC];
iniGuessX = [iniGuessEtaHat];



% - - - - - - - - - - - - - - - - - - - 
% Use the internal Matlab "Newton's Method
% - - - - - - - - - - - - - - - - - - - 
options = optimset('TolFun',1e-8,'TolX',1e-8,'jacobian','off','display','iter');
X = fsolve(@(X) solveFunct(X,x,N), iniGuessX,options);


% - - - - - - - - - - - - - - - - - - - 
% Retreive and plot results
% - - - - - - - - - - - - - - - - - - - 
etaHat = X(1:end-1);
eta = invHat(etaHat,x,N);
etax = d(eta,x,N);
c = X(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - - - - - - - - - - - - - - - - - - - 
% Compute the Fourier Coefficients
% - - - - - - - - - - - - - - - - - - - 
function uHat = hat(u,x,N)
dx = x(2)-x(1);
L = x(end)-x(1);
intTerm = exp(-1i*N*x).*u;      % create a matrix mesh of the Fourier integrand
uHat = sum(intTerm(:,2:end),2)*dx/L';   % This is the trapz. rule for periodic functions
end

% - - - - - - - - - - - - - - - - - - - 
% Compute the Fourier Series from Fourier Coefficients
% - - - - - - - - - - - - - - - - - - - 
function u = invHat(uHat,x,N)
sumTerm = exp(1i*N*x).*uHat;
u = sum(sumTerm,1);
if norm(imag(u))<1e-10
    u = real(u);
else
    disp('Complex function')
end
end


% - - - - - - - - - - - - - - - - - - - 
% Compute Derivative Spectrally
% - - - - - - - - - - - - - - - - - - - 
function uX = d(u,x,N)
uHat = hat(u,x,N);
uXHat = uHat.*(1i*N);
uX = invHat(uXHat,x,N);
end

% - - - - - - - - - - - - - - - - - - - 
% Compute the Nonlinear Eqns & Jacobian
% - - - - - - - - - - - - - - - - - - - 
function [F dF] = solveFunct(X,x,N);
Nn = max(N);
etaHat = X(1:end);
g = 1;
h = 0.5;
c = sqrt(g*tanh(h));
rho = 1;
P = surfaceToPressure;


eta = invHat(etaHat,x,N);
etaX = d(eta,x,N);

F = zeros(2*Nn+1,1);
F(Nn+1) = etaHat(Nn+1);

for j=1:2*Nn+1
    int = @(x)cos(j*x)*(c-sqrt(c^2-2.*(P-rho*g*h)));
    coeff = integral(int,0,2*pi,'ArrayValued',true);
    sigma(j,:) =(1/pi).*coeff.*cos(j*x).*cosh(j*(eta+h));
end

q = sum(sigma) - c + sqrt((c^2-2*rho*g*eta)./(1+(etaX).^2));
qHat = hat(q,x,N);
% qHat(end+1:numel(F))=0;
F = qHat;
% F = F + qHat;
% F(2*Nn+2) = sum(etaHat)-.0001;
%Since we don't know the maximum wave 
%amplitude, we'll end up running a feedback loop and adjusting this
 %accordingly. For the meantime, we'll leave it as it is
dF = jacobian(X,x,N);
end

% - - - - - - - - - - - - - - - - - - - 
% Calculate the Jacobian
% - - - - - - - - - - - - - - - - - - - 
function dF = jacobian(X,x,N)

    Nn = max(N);
    dF = zeros(2*Nn+1,2*Nn+1);

    etaHat = X(1:end);
    g = 1;
    h = 0.5;
    c=sqrt(g*tanh(h));
    rho = 1;
    P = surfaceToPressure;

    eta = invHat(etaHat,x,N);
    etaX = d(eta,x,N);

    for k = 1:1+2*Nn
        for j=1:2*Nn+1
            int = @(x)cos(j*x)*(c-sqrt(c^2-2.*(P-rho*g*h)));
            coeff = integral(int,0,2*pi,'ArrayValued',true);
            sigma(j,:) =(1/pi).*coeff.*cos(j*x).*sinh(j*(eta+h)).*(j*exp(1i*k*x));
        end

        dFdEta = sum(sigma) + (1/2)*((c^2-2*rho*g*eta)./(1+(etaX).^2)).^(-1/2)...
        *(((-2*rho*g*exp(1i*k*x)).*(1+(etaX).^2)-(c^2-2*rho*g*eta).*(2*etaX)...
        .*(i*k*exp(i*k*x)))/(1+(etaX).^2).^2);
    
        dF(k,:) = hat(dFdEta,x,N);
    end
end

