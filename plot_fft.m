function  plot_fft(sig, fs)

Y = fft(sig);
L = length(sig);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
fAxis = fs*(0:(L/2))/L;

plot(fAxis, P1)