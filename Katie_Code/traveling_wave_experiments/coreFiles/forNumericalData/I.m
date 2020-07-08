function fHat = I(eta,h);
M = length(eta);
N = (M-1)/2;

etaX = d(eta);
etaHat = hat(eta);

fHat = zeros(M,1);
for k = -N:N
    kInd = N+k+1;
    for n = -N:N
        nInd = n+N+1;
        if abs(k-n)<=N
            if (n==0)||(k==0)
                temp=0;
            else
                temp = hat((sinh(n*(eta + h))+i*etaX.*cosh(n*(eta + h))).^(-1),(k-n));
            end
            fHat(kInd) = fHat(kInd) - n*etaHat(nInd)*temp;
        end
        
    end
end