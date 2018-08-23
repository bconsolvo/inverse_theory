%% 1.1 
% Define and plot time series
ak=[20,-1,6,7,1,-5]; % Define time series values
t=[0:1:5]; % Define time variable
figure;
stem(t,ak); % vertical line plot
grid on;
hold on;
plot(t,ak); % continuous line plot
title('Time Series ak (causal)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-6 6]);

% shift the time series to the left
ak=[20,-1,6,7,1,-5]; % Define time series values
t=[-1:1:4]; % Define time variable
figure;
stem(t,ak); % vertical line plot
grid on;
hold on;
plot(t,ak); % continuous line plot
title('Time Series ak shifted (non-causal)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-6 6]);
