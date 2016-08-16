
clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC's and simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_t = 0.1;
t_start = 0;
t_end = 10;
t = t_start: delta_t: t_end;

x_IC = [1 1 1];

% Pre-allocate
x_OL = zeros(length(t), 3);
x_OL(1,:) = x_IC;
y_OL = zeros(length(t),1 );
y_OL(1) = x_IC(2);      % y = x2


for epoch = 2: length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Put the system in normal form, i.e. calculate xi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate u_eq with the switched Lyapunov algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate dV1_dot_du
    
    % Calculate dV2_dot_du
    
    % Compare dV1_dot_du and dV2_dot_du to choose the CLF
    
    % Calculate u_eq with the CLF of choice
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate u_s with the SMC algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate omega (the surface)
    
    % Calculate u_s
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply (u = u_eq + u_s) to the system and simulate for one time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate the open-loop system
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xy = simulate_sys( x_OL(epoch-1,:), y_OL(epoch-1), 0, delta_t);
    x_OL(epoch,:) = xy(1:3);    % First 3 elements are x
    y_OL(epoch) = xy(end);         % Final element is y
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
plot(t, x_OL)
legend('x_1','x_2','x_3','location','NorthWest')
xlabel('Time [s]')
ylabel('x')
title('x: Open Loop')

subplot(2,2,2)
plot(t, y_OL)
xlabel('Time [s]')
ylabel('x')
title('y: Open Loop')