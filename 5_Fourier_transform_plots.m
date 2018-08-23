%% Because ak is causal, its Fourier Transform should have real and imag parts that are Hilbert Transform pairs.
close all; clc; clear all;
ak=[20,-1,6,7,1,-5]; % Define time series values
t=[0:1:5]; % Define time variable
FTak=ifft(ak); % Taking the Fourier Transform of the time series ak. The IFFT "inverse fast fourier transform" is taken in our case because the convention in MATLAB is to put the normalizing factor in front of this function.
figure;
plot(t,real(FTak)); % Plotting the real part of the Fourier transform
hold on;
plot(t,imag(FTak)); % Plotting the imaginary part of the FT
legend('Real part of FT(a_k)','Imag part of FT(a_k)');
title('Real and Imaginary part of FT of a causal filter A_k');
xlabel('Frequency Bin');
ylabel('Amplitude');
%% Calculating the log amplitude spectrum of Ak and then the Hilbert Transform of it to get the phase
a1=abs(FTak); % log amplitude spectrum
h1=hilbert(a1); % calculating the analytic function from the log amplitude spectrum
p1=imag(h1); % The associated phase of Ak by taking the imaginary part of the analytic signal
figure;
plot(t,a1); % plotting the log of the amplitude spectrum
title('Log Amplitude Spectrum of A_k');
xlabel('Frequency Bin');
ylabel('Log Amplitude');

figure;
plot(t,p1); % plotting the Hilbert transform of the log of the amplitude spectrum (the phase)
title('Phase Spectrum of A_k by HT calculation');
xlabel('Frequency Bin')
ylabel('Phase');



