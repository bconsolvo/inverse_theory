%% Question 1.a)
t1=[0; 3.0]; % defining time vector
d1=[2;46]; % defining distance vector

% least squares solution to solve for velocity
v1=(inv(t1'*t1))*t1'*d1

% Let us now calculate the covariance of the model
Cd1=[5^2 0; 0 5^2] % covariance of data
Cm1=inv(t1'*inv(Cd1)*t1) %covariance of model

% Let us now calculate the root mean squared of the prediction errors
p1=d1-t1*v1 % prediction error vector
rms1=sqrt((2^2+0^2)/2) % rms of prediction errors


