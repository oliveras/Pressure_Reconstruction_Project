function ans = integrate(x,u);

x = [x; pi*2];
y = [u; u(1)];

ans = x(2)/2*(sum(y(1:end-1))+sum(y(2:end)));


