%%1.2 drawing the unit circle in the Argand plane
angle=linspace(0,2*pi,360); % Linearly spaced vector generates 360 points between 0 and 2 pi
x=cos(angle); % define x
y=sin(angle); % define y

% defining the roots
zx=[1.7611 -1.16212 -1.16212 0.381574 0.381574] % Real component of root
zy=[0 -0.787121 0.787121 -1.00364 1.00364] % Imaginary component of root


figure;
plot(x,y); % plot the circle
axis('equal'); % making the axes equal size so that the circle really looks like a circle
grid on;
set(gca,'XAxisLocation','origin','YAxisLocation','origin'); %moving axes to the centre
text(-0.5,0.1,'maximum delay'); % displaying where maximum delay occurs
text(1.2,1.2,'minimum delay'); % displaying where minimum delay occurs
xlabel('Re');
title('Minimum and Maximum Delay - Argand Plane Unit Circle');
ylabel('Im');
xlim([-2 2]); % setting limits for x-axis
ylim([-2 2]); % setting limits for y-axis
hold on; 
plot(zx,zy,'*');

