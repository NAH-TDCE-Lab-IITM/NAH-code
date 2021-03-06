function [P1]=fourier(sampling_rate,S)

L=length(S);
Y = fft(S);
P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
f = sampling_rate*(0:(L/2))/L;
% 
 figure(2)
 plot(f,P1) 
 title('Single-Sided Amplitude Spectrum of S(t)')
 xlabel('f (Hz)')
 ylabel('|P1(f)|')

end