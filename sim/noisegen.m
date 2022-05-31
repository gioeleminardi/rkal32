close all
clear all

RocketSimulation = import_data("/home/varkrid/code/cansat/rkal32/sim/RocketSimulation.csv", [2, Inf]);
RocketSimulation = rot90(RocketSimulation);

dt = (10^-2); % Time steps of 10ms  

t_length = 18.14;

t = [0:dt:t_length]; % Create Time vector

%% Measurement noise
% Position corrupted with Gaussian noise with standard deviation = 5.985
s_n = RocketSimulation(1,:)+(5.985.*randn(1,length(t)));
% Velocity corrupted with Gaussian noise with standard deviation = 2
v_n = RocketSimulation(3,:)+(2.*randn(1,length(t)));
% Acceleration corrupted with Gaussian noise with standard deviation = 0.0346
a_n = RocketSimulation(2,:)+(0.0346.*randn(1,length(t)));

figure;
plot(t,s_n) % Plot noisy position
title('Noisy Position')

figure;
plot(t,v_n) % Plot noisy velocity
title('Noisy Velocity')

figure;
plot(t,a_n) % Plot noisy acceleration
title('Noisy Acceleration')

out_matrix(:,1) = t;
out_matrix(:,2) = v_n;
out_matrix(:,3) = a_n;
out_matrix(:,4) = s_n;

writematrix(out_matrix, 'out_matrix.csv', 'Delimiter',',')
