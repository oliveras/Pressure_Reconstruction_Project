function [p Phi PhiX] = generatePressureData(eta,c,h,g)

M = length(eta); N = (M-1)/2; x = (0:2*N)'/M*2*pi;
kv = (-N:N)'; etaX = d(eta);

b = hat(sqrt(c^2-2*g*eta)./sqrt(1 + etaX.^2));

for n=-N:N
    rowInd = n + N + 1;
    for k=-N:N
        colInd = k + N + 1;
        if abs(n-k)<=N/2
            A(rowInd,colInd) =  hat(cosh(k*eta+k*h),n-k);
        end
    end
end

PhiHat = A\b; 

Phi = invHat(PhiHat);  PhiX = (Phi);

p = .5*(c^2 - (PhiX).^2); 