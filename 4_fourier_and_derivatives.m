
% Exercises in Fourier Transforms 2
% Completed during MS Geophysics at the University of Western Ontario
% MATLAB code by Benjamin Consolvo and Gerhard Pratt
% Updated in 11/2015

clear all;
clc;
close all;

%%  Fourier transform of the Boxcar function

M = 5;
box = zeros(100,1);
box(1:M)=1;
N = length(box);            % Find the length of time series
X = zeros(N,1);           % Set up an N-by-1 matrix of zeros

X=ifft(box); % inverse Fourier transform

int=0:N;
freq=-max(int/2):1:max(int/2);
freq=freq';

g1=zeros(N+1,1);

g1(1:50)=X(51:100);
g1(51:100)=X(1:50);
g1(101)=0.0099;

ampspect=abs(g1);

% Plot 4a)
ps=angle(g1);
figure;
plot(freq,ps,'linewidth',2);
title('Phase Spectrum of An');
xlabel('Frequency (Hz)');
ylabel('Phase Spectrum');
legend('M=5');

% Plot 4b)
figure;
plot(freq,ampspect,'linewidth',2);                
title('Amplitude Spectrum of An');
xlabel('Frequency (Hz)');
ylabel('Amplitude Spectrum');
legend('M=5');

%% Testing different values of M
M = 5;
box = zeros(100,1);
box(1:M)=1;
N = length(box);            % Find the length of time series
X = zeros(N,1);           % Set up an N-by-1 matrix of zeros

X=ifft(box);

int=0:N;
freq=-max(int/2):1:max(int/2);
freq=freq';

g1=zeros(N+1,1);

g1(1:50)=X(51:100);
g1(51:100)=X(1:50);
g1(101)=0.0099;

ampspect=abs(g1);
% Plot 4c)
figure;
plot(freq,ampspect,'linewidth',2);                
title('Amplitude Spectrum of An');
xlabel('Frequency (Hz)');
ylabel('Amplitude Spectrum');
legend('M=5');
hold on;

M = 10;
box = zeros(100,1);
box(1:M)=1;
N = length(box);            % Find the length of time series
X = zeros(N,1);           % Set up an N-by-1 matrix of zeros

X=ifft(box);

int=0:N;
freq=-max(int/2):1:max(int/2);
freq=freq';

g1=zeros(N+1,1);

g1(1:50)=X(51:100);
g1(51:100)=X(1:50);

ampspect=abs(g1);

plot(freq,ampspect,'r','linewidth',2);                
title('Amplitude Spectrum of An');
xlabel('Frequency (Hz)');
ylabel('Amplitude Spectrum');
hold on;

M = 20;
box = zeros(100,1);
box(1:M)=1;
N = length(box);            % Find the length of time series
X = zeros(N,1);           % Set up an N-by-1 matrix of zeros

X=ifft(box);

int=0:N;
freq=-max(int/2):1:max(int/2);
freq=freq';

g1=zeros(N+1,1);
g1(1:50)=X(51:100);
g1(51:100)=X(1:50);

ampspect=abs(g1);

plot(freq,ampspect,'g','linewidth',2);                
title('Amplitude Spectrum of An');
xlabel('Frequency (Hz)');
ylabel('Amplitude Spectrum');
legend('M=5; N/M=20','M=10; N/M=10','M=20; N/M=5');
set(gca,'xtick',[-100:5:100]);

%% Work with Gaussian curve
% 4d) Gaussian in time domain
t=[-5:.1:5];
alpha=1;
gt=exp(-alpha*(t.^2));
figure;
plot(t,gt);
title('Gaussian Curve in Time');
xlabel('Time (s)');
ylabel('Amplitude');
%eb=0.707*ones(size(gt));
%errorbar(gt,eb);

% 4e) Gaussian in frequency domain
fs=1/0.1;
N=length(gt);
n=[-50:1:50];
omega=n*2*pi/(N*0.1);
%gf=fft(gt);
f=(n*fs/N);
gf=(sqrt(pi/alpha))*exp(-(omega.^2)/(4*alpha));
figure;
plot(f,gf);
title('Gaussian Curve in Frequency');
xlabel('Frequency Bin');
ylabel('Amplitude');


%%  Derivative theorem test
close all;
clear all;
clc;
x=[0:(pi/16):8*pi]; 
an=cos(x); % Defining time series

% An=ifft(an); % Takes the forward DFT in Gerhard's notation
% First, let us take the FT of an.
N = length(an);            % Find the length of time series
An = zeros(N,1);           % Set up an N-by-1 matrix of zeros

for n = 1:N
    for k = 1:N
        compl=i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        An(n)=An(n)+(1/N)*an(k)*exp(compl); %Formula for DFT
    end
end

% Code to carry out the derivative of An

N=length(An);
anderiv = zeros(N,1);
for n = 1:N
    for k = 1:N
        omegan=(2*pi*(n-1))/(N*(pi/16));
        compl=-i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        anderiv(k)=anderiv(k)+0.066*i*omegan*An(n)*exp(compl); %Formula for inverse DFT to get ak back
        % Had to multiply by 0.066 to get my derivative to match.
    end
end

% Plot 4f)
figure;
plot(x,real(anderiv),'r',x,-sin(x),'g');
title('Derivative of Time Series ak');
xlabel('Time (s)');
ylabel('Amplitude');
hold on;

% ii) approximate derivative of ak - "first difference"

t=pi/16;
akderiv=zeros(length(an),1);
akderiv=akderiv';

for i=1:length(an)
    if i+1>length(an)
        an(i+1)=1;
    end
    akderiv(i)=akderiv(i)+((an(i+1)-an(i))/t);
    
end
x2=x(1:129);
akderiv=akderiv(1:129);

plot(x2,akderiv,'b')
title('Derivative of Time Series ak');
xlabel('Time (s)');
ylabel('Amplitude');
legend('0.066*akderiv','-sin(x)','first difference derivative')





%% Demonstrating the Poisson relationship: horizontal position

close all;
clear all;
clc;
% Assume a value of unity for constants
load data/gprofile.dat;
grav=gprofile';
x=[0:1.4:336];
an=grav; % Defining time series

% An=ifft(an); % Takes the forward DFT
% First, let us take the FT of an.
N = length(an);            % Find the length of time series
An = zeros(N,1);           % Set up an N-by-1 matrix of zeros

for n = 1:N
    for k = 1:N
        compl=i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        An(n)=An(n)+(1/N)*an(k)*exp(compl); %Formula for DFT
    end
end

% Code to carry out the derivative of an

N=length(An);
anderiv = zeros(N,1);
for n = 1:N
    for k = 1:N
        omegan=(2*pi*(n-1))/(N*(pi/16));
        compl=-i*2*pi*(n-1)*(k-1)/N; % Calculating the complex exponent
        anderiv(k)=anderiv(k)+0.066*i*omegan*An(n)*exp(compl); %Formula for inverse DFT to get ak back
        % I had to multiply by 0.066 to get my derivative to match.
    end
end


figure;
plot(x,grav,'b',x,real(anderiv),'r');
title('Gravity Field and Pseudo Magnetic Field');
xlabel('Distance (km)');
ylabel('Potential Field Amplitude');
legend('Gravity Profile','Dgrav = Pseudo Magnetic')
hold on;
