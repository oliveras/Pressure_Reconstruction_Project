function f = zeroMode(f)

fHat = fft(f);
fHat(1) = 0;
f = ifft(fHat);
