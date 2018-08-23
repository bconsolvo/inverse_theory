%%1.2 drawing a circle for the minimum phase diagram
clc; clear all; close all; % clearing command line, Clearing variables, closing figure windows
angle=linspace(0,2*pi,360); % Linearly spaced vector generates 360 points between 0 and 2 pi
x=2*cos(angle)+2.5; % define x
y=2*sin(angle); % define y


figure;
plot(x,y); % plot the circle
axis('equal'); % making the axes equal size so that the circle really looks like a circle
xlabel('Real');
ylabel('Imaginary');
set(gca,'XAxisLocation','origin','YAxisLocation','origin'); %moving axes to the centre
set(gca,'YTick',[]); % Turning numbers on grid off
set(gca,'XTick',[]); % Turning numbers on grid off
title('Minimum Phase');
xlim([-2 8]); % setting limits for x-axis
ylim([-3 4]); % setting limits for y-axis
