function uFilter = filterData(u,nn,mm)

if nargin>2
    M = mm/2;
    N = nn/2;
    uFilter = interpft(u,mm);
    uHat = fftshift(fft(uFilter));
    tails = [ones(N-1,1); sin(pi/2*(((N:M)'-N)/(M-N)+1))];
    vals = [flipud(tails);1;tails(1:end-1)];
    
    uHat = uHat.*vals; 
    uFilter = ifft(ifftshift(uHat));
    
else
    uFilter = interpft(u,nn);
end

%uFilter = interpft(uFilter,length(u));