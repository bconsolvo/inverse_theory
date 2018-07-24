% Autocorrelation function
% Completed during MS Geophysics at the University of Western Ontario
% MATLAB code by Benjamin Consolvo
% Updated in 12/2015

% Autocorrelation of a function F(z)
F=[-256 1536 -3840 5376 -4704 2688 -1008 240 -33 2]; %Defining time series F
x=[1 2 3 4 5 6 7 8 9 10] % Defining position in time series for plot
Finv = fliplr(F); %Calculating the time reverse of F
ACF=conv(F,Finv); %Autocorrelation of F
disp(ACF) %Displays autocorrelation

%% Time Series with the same autocorrelation
G=[512 -2304 4608 -5376 4032 -2016 672 -144 18 -1]; %Defining time series G
xone=[1 2 3 4 5 6 7 8 9 10] % Defining position in time series for plot
Ginv=fliplr(G); %Calculating the time reverse of G
ACG=conv(G,Ginv); %Autocorrelation of G
disp(ACG) %Displays autocorrelation

% plotting both time series dipole factors (or polynomial coefficients)
figure;
close all;
plot(x,F,xone,G);
xlabel('Time Series Position');
ylabel('Amplitude (or dipole factors)');
title('Dipole Factors of two Time Series with the same autocorrelation');
s1='F(z) time series';
s2='G(z) time series';
legend(s1,s2);
