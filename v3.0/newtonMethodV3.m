function returnVal = newtonMethodV3(iniGuessEtaHat,iniGuessC, a, h, g, rho, N, L, x)

    iniGuessX = [iniGuessEtaHat;iniGuessC];
    % - - - - - - - - - - - - - - - - - - - 
    % Use the internal Matlab "Newton's Method"
    % - - - - - - - - - - - - - - - - - - - 
    options = optimset('TolFun',1e-6,'TolX',1e-6,'jacobian','on','display','iter','MaxIter',50000,'MaxFunEvals',50000);
    [X,fval,exitflag,output,Jacobian] = fsolve(@(X) solveFunct(X,x,N), iniGuessX,options);

    etaHat = X(1:end-1);
    eta = invHat(etaHat,x,N);
    c = X(end);

    returnVal = {eta,c,Jacobian};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function uHat = hat(u,x,N)
        dx = x(2)-x(1);
        L = x(end)-x(1);
        intTerm = cos(N*x).*u;      % create a matrix mesh of the Fourier integrand
        uHat = 2*sum(intTerm(:,2:end),2)*dx/L';   % This is the trapz. rule for periodic functions
    end

    % - - - - - - - - - - - - - - - - - - - 
    % Compute the Fourier Series from Fourier Coefficients
    % - - - - - - - - - - - - - - - - - - - 
    function u = invHat(uHat,x,N)
        L = x(end)-x(1);
        uHat(1) = uHat(1)/2;
        sumTerm = cos(N*x).*uHat;
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

        F = zeros(Nn+2,1);
        for n = 0:Nn
            if n==0
                F(1) = 0;
            else
                F(n+1) = hat((g/n)*cosh(n*(eta+h))+(c^2-2*g*eta).*sinh(n*(eta+h)),x,n);
            end
        end

        F(end) = sum(etaHat)-a;

        dF = jacobian(X,x,N);
    end

    % - - - - - - - - - - - - - - - - - - - 
    % Calculate the Jacobian
    % - - - - - - - - - - - - - - - - - - - 
    function dF = jacobian(X,x,N)

        Nn = max(N);
        dF = zeros(Nn+2,Nn+2);

        etaHat = X(1:end-1);
        c = X(end);
        g = 1;
        h = 0.5;
        L = x(end)-x(1);

        eta = invHat(etaHat,x,N);
 
        for n = 0:Nn
            if n==0
                dF(1,1) = 1;
            else
                for k = 0:Nn
                    dFdEta = cos(k*x).*(n*(c^2-2*g*eta).*cosh(n*(eta+h))-g*sinh(n*(eta+h)));
                    dF(n+1,k+1) = hat(dFdEta,x,n)*(abs(n-k) <= Nn);
                end
            end
            dF(n+1,Nn+2) = hat(2*c.*sinh(n*(eta+h)),x,n);
        end
        dF(Nn+2,1:Nn+1) = ones(1,Nn+1);
        dF(Nn+2, Nn+2) = 0;
    end

end