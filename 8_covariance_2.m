%% Question 1.b)
clear all;
clc; close all;

% i) Calculating for model where d=vt
t2=[0; 1.5; 3.0]; % defining time vector
d2=[2; 12; 46]; % defining distance (data) vector
v2=(inv(t2'*t2))*t2'*d2
Cd2=[5^2 0 0; 0 5^2 0; 0 0 5^2] % covariance of data
Cm2=inv(t2'*inv(Cd2)*t2) %covariance of model
sigma2=sqrt(Cm2) % standard deviation of model
p2=d2-t2*v2 % prediction error vector
rms2=sqrt((2^2+(-8.8)^2+(4.4)^2)/3) % rms of prediction errors



% ii) Calculations for model where d=vt+.5at^2

ta=[0 .5*(t2(1))^2; 1.5 .5*(t2(2))^2; 3.0 .5*(t2(3))^2]   % defining time and acceleration matrix
v22=(inv(ta'*ta))*ta'*d2 % least squares solution to solve for velocity

% Let us now calculate the covariance of the model
Cd22=[5^2 0 0; 0 5^2 0; 0 0 5^2] % covariance of data
Cm22=inv(ta'*inv(Cd22)*ta) %covariance of model

% Standard deviation of velocity
sigma22=[sqrt(Cm22(1)); sqrt(Cm22(4))]

% standard deviation calculation that did not work :(
% sd2=5/(sqrt(0^2+1.5^2+3^2))

p22=d2-ta*v22 % prediction error vector
rms22=sqrt((2^2+0^2+0^2)/3) % rms of prediction errors