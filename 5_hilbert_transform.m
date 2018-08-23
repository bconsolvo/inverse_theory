%% Hilbert Transform 
% Example of HT with a sin function
close all;
clc;
clear all;
x=0:(pi/30):2*pi % x-interval defined
y=sin(x); % Defining the series
plot(x,y);
ht=hilbert(y); % Taking the Hilbert Transform of ak - this creates the analytical signal with a real and imaginary part
hold on;
plot(x,imag(ht));
legend('sin(x)','Imag part (HT)')
title('Hilbert Transform of a sin function')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%% Example of HT with constructed causal, minimum phase function (ak)
close all; clc; clear all;
% Define and plot time series
ak=[20,-1,6,7,1,-5]; % Define time series values
t=[0:1:5]; % Define time variable
figure;
%stem(t,ak); % vertical line plot - not needed here
grid on;

plot(t,ak); % continuous line plot
title('Time Series ak (minimum phase and causal)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-2 8]);
hold on;
Hak=hilbert(ak); % Taking the Hilbert Transform of ak - this creates the analytical signal with a real and imaginary part
hold on;
plot(t,imag(Hak));
legend('time series a_k','HT pair of a_k')
title('Hilbert Transform of a time series a_k')
xlabel('Time (s)')
ylabel('Amplitude')


%% Example with real data (gprofile.dat)

close all;
clear all;
clc;
load gprofile.dat; % Loading gravity data
grav=gprofile';
x=[0:1.4:336];
an=grav; % Defining time series for gravity

HTg=hilbert(an);

plot(x,an);
hold on;
plot(x,imag(HTg));
title('Hilbert Transform of Gravity Profile Data');
xlabel('Distance (km)');
ylabel('Potential Field Amplitude');
legend('Gravity Profile','Hilbert Transform (imaginary part)')



