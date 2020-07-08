function output = solveFuncNew(X,x,g,amp,sigma,h)

eta = X(1:end-1); c = X(end);N = (length(eta)-1)/2;
[cc cV] = getC(eta,h,g);
firstEqn = hat(c^2 - cV.^2);
firstEqn(N+1)=0;

finalEqn = abs(amp - .5*(max(eta) - min(eta)));

output = [firstEqn; finalEqn];
