clear all
close all
clc


N = 32; M = 2*N+1;
x = (0:2*N)'/M*2*pi;

p.N = N; p.M = M; p.x = x;



lambda = 0.001; p.lambda = lambda;
eta0 = lambda*cos(x);
etaHat0 = hat(eta0);
c0 = sqrt(tanh(1));
X0 = [etaHat0; c0];


options = optimset(...
    'display',      'off',...
    'tolfun',       1e-15,...
    'tolx',         1e-15,...
    'jacobian',     'on');% , ...

Flambda = zeros(M+1,1);
Flambda(end) = -1;
lambdaVals = [lambda];
tic
X = fsolve(@(X) solveFunc(X,p), X0,options);
if norm(imag(X))<1e-15
    X = real(X);
end
subplot(3,1,1)
plot(x,real(invHat(X(1:M))))
axis tight
subplot(3,1,2)
stem(-N:N,X(1:M))


subplot(3,1,3)
plot([c0 X(end)],[0,lambda],'*','MarkerFaceColor','k','color','k')
hold on

dLambda = .005;
for contIter = 2:50
    [F DF] = solveFunc(X,p);
    lambda = lambda + dLambda; p.lambda = lambda;
    lambdaVals(contIter) = lambda;
    X0 = X - (dLambda)*(DF\Flambda);
    
    X = fsolve(@(X) solveFunc(X,p), X0,options);
    if(norm(imag(X))<1e-15)
        X = real(X);
    end
    etaHat = X(1:M);
    eta = real(invHat(etaHat));
    
    subplot(3,1,1)
    plot(x,eta);
    axis tight
    subplot(3,1,2)
    stem(-N:N,X(1:M))

    subplot(3,1,3)
    plot(X(end),lambda,'*','MarkerFaceColor','k','color','k')
    hold on
    pause(.1)
end
toc

function [F DF] = solveFunc(X,p)

    N = p.N; M = p.M; x = p.x; lambda = p.lambda;
    F = zeros(M+1,1); 

    etaHat = X(1:M);
    c = X(end);
    eta = invHat(etaHat);

    kV = (-N:N)';
    T = tanh(kV)./kV; T(N+1) = 0;

    F(1:M) = hat(sqrt(c^2 - 2*eta)).*T + hat(eta.*sqrt(c^2 - 2*eta));
    F(N+1) = etaHat(N+1);
    F(end) = sum(etaHat) - lambda;

    F = real(F);
    DF = zeros(M+1,M+1);

    for j=1:M
        if j==N+1
            DF(j,j) = 1/(2*pi);
        else
            k = j - N - 1;
            term1 = -hat(exp(1i*k*x)./sqrt(c^2 - 2*eta)).*tanh(k)/k;
            term2 = -hat(exp(1i*k*x)./sqrt(c^2 - 2*eta).*eta);
            term3 = hat(exp(1i*k*x).*sqrt(c^2 - 2*eta));
            DF(j,1:M) = 1/(2*pi)*(term1 + term2 + term3)';
        end
    end
    DF(1:M,M+1) = 1/(2*pi)*(hat(c./sqrt(c^2-2*eta)).*T + hat(c./sqrt(c^2-2*eta).*eta));
    DF(N+1,M+1) = 0;
    DF(M+1,1:M) = 1;
    DF(M+1,M+1) = 0;
    DF = real(DF);
end



function uHat = hat(u);
    uHat = fftshift(fft(u))/length(u);
end

function u = invHat(uHat);
    u = real(ifft(ifftshift(uHat))*length(uHat));
end