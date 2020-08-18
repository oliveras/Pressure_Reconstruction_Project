clear all
close all
clc

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
Nn = 32;                        % number of Fourier Modes
N = (-Nn:Nn)';                  % a vector of modes

% - - - - - - - - - - - - - - - - - - - 
% Initial Guess Data
% - - - - - - - - - - - - - - - - - - - 
a = 0.01:0.01:0.5; % use increments of 0.01 for a values
X = cell(size(a,2)+1,2); % store past iterations for initial guesses

iniGuessEta = a(1)*cos(x);
iniGuessEtaHat = hat(iniGuessEta,x,N);
iniGuessC = sqrt(g*tanh(h));
% iniGuessX = [iniGuessEtaHat;iniGuessC];

X{1,1} = iniGuessEtaHat;
X{1,2} = iniGuessC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(a)
    [eta c] = newtonMethod(X{i,1},X{i,2}, a(i), h, g, rho, N, L, x);
    X{i+1,1} = hat(eta,x,N);
    X{i+1,2} = c;
    axes_options = {'interpreter','latex','fontsize',12};
    subplot(3,1,1)
    stem(N, real(hat(X{i,1},x,N)))
    xlabel('$N$',axes_options{:})
    title('Fourier Coefficients of $\eta$',axes_options{:})

    subplot(3,1,2)    
    if i > 1
        plot(x,invHat(X{i+1,1},x,N), '-.',x,invHat(X{i,1},x,N))
    else
        plot(x,invHat(X{i+1,1},x,N))
    end
    
    subplot(3,1,3)
    scatter([X{2:end,2}],a(1:length([X{2:end,2}])), 'filled')
    
    pause(2.0)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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