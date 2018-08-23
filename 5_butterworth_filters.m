%% 3. b)
% Plotting poles for n = 2,4,6
clear all;
n=2;
m=[-n:1:n];
omegac=1;
for k=1:length(m);
    u(k) = omegac*exp((i*pi*(2*m(k)+1+n)/(2*n)));
end 


figure;
plot(real(u),imag(u),'*');
set(gca,'XAxisLocation','origin','YAxisLocation','origin');
title('Butterworth Filter Poles n=2, \omega_c=1');
xlabel('Real');
ylabel('Imaginary');
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
grid on;
hold on;
plot(real(u),imag(u));

clear all;
n=4;
m=[-n:1:n];
omegac=1;
for k=1:length(m);
    u(k) = omegac*exp((i*pi*(2*m(k)+1+n)/(2*n)));
end 


figure;
plot(real(u),imag(u),'*');
set(gca,'XAxisLocation','origin','YAxisLocation','origin');
title('Butterworth Filter Poles n=4, \omega_c=1');
xlabel('Real');
ylabel('Imaginary');
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
grid on;
hold on;
plot(real(u),imag(u));


clear all;
n=6;
m=[-n:1:n];
omegac=1;
for k=1:length(m);
    u(k) = omegac*exp((i*pi*(2*m(k)+1+n)/(2*n)));
end 


figure;
plot(real(u),imag(u),'*');
set(gca,'XAxisLocation','origin','YAxisLocation','origin');
title('Butterworth Filter Poles n=6, \omega_c=1');
xlabel('Real');
ylabel('Imaginary');
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
grid on;
hold on;
plot(real(u),imag(u));
