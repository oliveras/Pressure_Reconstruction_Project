function [eta c] = newtonMethod(iniGuessEtaHat,iniGuessC, a, h, g, rho, N, L, x)

    iniGuessX = [iniGuessEtaHat;iniGuessC];
    % - - - - - - - - - - - - - - - - - - - 
    % Use the internal Matlab "Newton's Method"
    % - - - - - - - - - - - - - - - - - - - 
    options = optimset('TolFun',1e-6,'TolX',1e-6,'jacobian','on','display','iter');
    X = fsolve(@(X) solveFunct(X,x,N), iniGuessX,options);

    % - - - - - - - - - - - - - - - - - - - 
    % Retreive and plot results
    % - - - - - - - - - - - - - - - - - - - 
    etaHat = X(1:end-1);
    eta = invHat(etaHat,x,N);
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
    % Compute the Nonlinear Eqns & Jacobian
    % - - - - - - - - - - - - - - - - - - - 
    function [F dF] = solveFunct(X,x,N);
        Nn = max(N);
        etaHat = X(1:end-1);
        c = X(end);
        g = 1;
        h = 0.5;

        eta = invHat(etaHat,x,N);

        F = zeros(2*Nn+2,1);
        for n = -Nn:Nn

            if n==0
                F(Nn+1) = etaHat(Nn+1);
            else
                F(Nn+n+1) = hat((g/n)*cosh(n*(eta+h))+(c^2-2*g*eta).*sinh(n*(eta+h)),x,n);
            end
        end

        F(2*Nn+2) = sum(etaHat)-a;

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

        for n = -Nn:Nn
            if n==0
                dF(Nn+1,Nn+1) = 1;
            else
                for k = -Nn:Nn
                    dFdEta = exp(1i*k*x).*(n*(c^2-2*g*eta).*cosh(n*(eta+h))-g*sinh(n*(eta+h)));
                    dF(n+Nn+1,k+Nn+1) = hat(dFdEta,x,n)*(abs(n-k) <= Nn);
                end
            end
            dF(n+Nn+1,2*Nn+2) = hat(2*c.*sinh(n*(eta+h)),x,n);
        end
        dF(2*Nn+2,1:2*Nn+1) = ones(1,2*Nn+1);
        dF(2*Nn+2, 2*Nn+2) = 0;
    end

end