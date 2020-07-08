function P = surfaceToPressure(~)
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
    iniGuessX = [iniGuessEtaHat;iniGuessC];


    % - - - - - - - - - - - - - - - - - - - 
    % Use the internal Matlab "Newton's Method
    % - - - - - - - - - - - - - - - - - - - 
    options = optimset('TolFun',1e-8,'TolX',1e-8,'jacobian','on','display','iter');
    X = fsolve(@(X) solveFunct(X,x,N), iniGuessX,options);


    % - - - - - - - - - - - - - - - - - - - 
    % Retreive and plot results
    % - - - - - - - - - - - - - - - - - - - 
    etaHat = X(1:end-1);
    eta = invHat(etaHat,x,N);
    etax = d(eta,x,N);
    c = X(end);


    % - - - - - - - - - - - - - - - - - - - 
    % Broke the qx term into two parts
    % - - - - - - - - - - - - - - - - - - - 
    sqrtTerm = sqrt((c^2 - 2*g.*eta).*(1 + etax.^2));
    qx = c + sqrtTerm; 


    % - - - - - - - - - - - - - - - - - - - 
    % Calculate the fk
    % - - - - - - - - - - - - - - - - - - - 
    fk = FC(etax, eta, qx, c, h, x, Nn);
    fk(Nn+1) = c - hat(sqrtTerm,x,0);  % modified the zeroth coefficient



    % - - - - - - - - - - - - - - - - - - - 
    % Added better plotting options
    % - - - - - - - - - - - - - - - - - - - 
    
    
%     subplot(2,1,1)
%     stem(N, real(fk))
%     xlabel('$N$',axes_options{:})
%     title('Fourier Coefficients of $\hat\phi_x(x,-h)$',axes_options{:})
% 
%     subplot(2,1,2)
    fkInv = invHat(fk, x, N);
    P = -rho*((fkInv -c).^2 - c.^2)/2;
%     plot(x, P,'-')
%     hold on
%     plot(x,rho*g*eta)
%     plot(x,eta)
%     leg = legend('$p- gh$','$\rho g \eta$','$\eta$','location','SouthEast');
%     set(leg,axes_options{:});
%     axis tight
%     xlabel('$x$',axes_options{:})


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % - - - - - - - - - - - - - - - - - - - 
    % Calculate the Fourier Coefficients of phi_x(x,-h)
    % - - - - - - - - - - - - - - - - - - - 
    function fk = FC(etax, eta, qx, c, h, x, Nn)
    fk = zeros(2*Nn+1, 1);
        for k = -Nn:Nn
            int = c.*etax.*sinh(k*(eta + h)) - 1i*qx.*cosh(k*(eta + h));
            fk(Nn+k+1) = -1i*hat(int, x, k);  % added a factor of 1i?  Is this correct?
        end
    end

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
    etaHat = X(1:end-1);
    c = X(end);
    g = 1;
    h = 0.5;



    eta = invHat(etaHat,x,N);
    etaX = d(eta,x,N);

    sqrtTerm = sqrt((c^2 - 2*g*eta).*(1 + etaX.^2));
    F = zeros(2*Nn+2,1);
    for n = -Nn:Nn

        if n==0
            F(Nn+1) = etaHat(Nn+1);
        else
            F(Nn+n+1) = hat(sqrtTerm.*sinh(n*(eta+h)),x,n);
        end
    end

    F(2*Nn+2) = sum(etaHat)-.0001;
    dF = jacobian(X,x,N);
    end

    % - - - - - - - - - - - - - - - - - - - 
    % Calculate the Jacobian
    % - - - - - - - - - - - - - - - - - - - 
    function dF = jacobian(X,x,N)

    Nn = max(N);
    dF = zeros(2*Nn+2,2*Nn+2);

    etaHat = X(1:end-1);
    c = X(end);
    g = 1;
    h = 0.5;

    eta = invHat(etaHat,x,N);
    etaX = d(eta,x,N);
    sqrtTerm = sqrt((c^2 - 2*g*eta).*(1 + etaX.^2));


    for n = -Nn:Nn
        if n==0
            dF(Nn+1,Nn+1) = 1;
        else
            for k = -Nn:Nn
                sinhTerm = sinh(n*(eta+h));
                dFdEta = (sinhTerm./sqrtTerm.*(-g.*(1+etaX.^2)+(c^2 - 2*g*eta).*etaX.*1i*k)...
                    + cosh(n*(eta+h)).*sqrtTerm.*n).*exp(1i*k*x);
                dF(n+Nn+1,k+Nn+1) = hat(dFdEta,x,n);
            end
        end
        dF(n+Nn+1,2*Nn+2) = hat(c./sqrtTerm.*(1+etaX.^2).*sinhTerm,x,n);
    end
    dF(2*Nn+2,1:2*Nn+1) = ones(1,2*Nn+1);

    end

end 