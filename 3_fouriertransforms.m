% Fourier Transforms, Amplitude and Phase Spectra
% Completed during MS Geophysics at the University of Western Ontario
% MATLAB code by Benjamin Consolvo
% Updated in 10/2015

close all;
clc;
clear all;
x=[0 1 2 3] % Frequency bins 


a1=complex(-.25,-.25)
a2=.5;
a3=complex(-.25,.25);
an=[0,a1,a2,a3]; % The complex function


ampspect=abs(an); %absolute value of the elements of an. When complex, 
                  % abs(an) is the complex modulus (magnitude) of
                  % the elements of an.
phasespect=angle(an); % Finding the phase spectrum of An

% a)
figure;
plot(x,phasespect);              % Plotting the amplitude spectrum of An
title('Phase Spectrum of An');
xlabel('Frequency Bin');
ylabel('Phase Spectrum');

% b)
figure;
plot(x,ampspect);                 % Plotting the phase spectrum of An
title('Amplitude Spectrum of An');
xlabel('Frequency Bin');
ylabel('Amplitude Spectrum');



%% Writing code to carry out the Discrete Fourier Transform (DFT)
a = [0 -1 1 0];           % Set time series to #1 question
N = length(a);            % Find the length of time series
X = zeros(N,1);           % Set up an N-by-1 matrix of zeros

for n = 1:N
    for k = 1:N
        compl=i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        X(n)=X(n)+(1/N)*a(k)*exp(compl); %Formula for DFT
    end
end
X;

% ii) 
% Writing code to carry out the inverse DFT

aa = zeros(N,1);
for n = 1:N
    for k = 1:N
        compl=-i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        aa(k)=aa(k)+X(n)*exp(compl); %Formula for inverse DFT to get ak back
    end
end
aa;

%% Taking a Fourier transform of the Boxcar function
% Writing code to carry out the DFT
box = zeros(100,1);
box(1:10)=1;
N = length(box);            % Find the length of time series
X = zeros(N,1);           % Set up an N-by-1 matrix of zeros

for n = 1:N
    for k = 1:N
        compl=i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        X(n)=X(n)+(1/N)*box(k)*exp(compl); %Formula for DFT
    end
end
X;
% Plotting the amplitude spectrum of a boxcar function An

f=zeros(100,1);
for q=1:99
    f(q+1)=f(q)+0.1;
end
 
f=f(1:51);
ampspect2=abs(X(1:51));

% c) 
figure;
plot(f,ampspect2);                
title('Amplitude Spectrum of An');
xlabel('Frequency (Hz)');
ylabel('Amplitude Spectrum');

                     
%%

% Read the time series "filteredtrace.dat" into MATLAB by using a dlmread
% command.
ex = dlmread('3_filteredtrace.dat');
t=zeros(1024,1);        % make a vector of length 1024 to put time values in
for y=1:1023
    t(y+1)=t(y)+0.002;  %convert to time in seconds 2 / 1000 = 0.002
end

% d)
figure;
plot(t,ex);                
title('Time Series filteredtrace.dat');
xlabel('Time (s)');
ylabel('Amplitude');

% ex is already our defined time series
N = length(ex);            % Find the length of time series
X = zeros(N,1);           % Set up an N-by-1 matrix of zeros

for n = 1:N
    for k = 1:N
        compl=i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        X(n)=X(n)+(1/N)*ex(k)*exp(compl); %Formula for DFT
    end
end
X;

%% Power spectrum is the square of the amplitude spectrum

t=0.002;
ff=zeros(N,1);
qu=[1:N];
qu=qu';
ff=qu/(N*t);

% e) Power spectrum from 0 to 500 Hz
pwrspect2=(abs(X)).^2;
figure;
plot(ff,pwrspect2);                
title('Power Spectrum of filteredtrace.dat');
xlabel('Frequency (Hz)');
ylabel('Amplitude Spectrum Squared');

% f) Power spectrum from 0 to250 Hz
figure;
plot(ff(1:N/2),pwrspect2(1:N/2));
title('Power Spectrum of filteredtrace.dat');
xlabel('Frequency (Hz)');
ylabel('Amplitude Spectrum Squared');



%% Apply phase shift in the frequency domain to advance the 
%seismic trace 5 ms in time, and then carry out the inverse DFT.  
%Produce plots to verify the time shift.


exft=fft(ex); % fourier transform of original time series
n=zeros(N,1); % creating a row vector of 0s
n=[1:1:512,-512:1:-1]; % This is the key part: we must index our n 
%appropriately to get a real valued time series when we take the IFFT.
%Here, we create a Hermite-symmetric complex exponential.
n=n'; % transpose of n

compl3=exp(complex(0,-2*pi*(n)*2.5/(N)));
ex_shifted=exft.*(compl3);

% g)
times=ifft(ex_shifted);
figure;
plot(t*1000,real(times));
hold on;
plot(t*1000,real(ex));
title('Original and Shifted Time Series filteredtrace.dat');
xlabel('Time (ms)');
ylabel('Amplitude');
legend('5ms shifted time series','Original Time Series');

