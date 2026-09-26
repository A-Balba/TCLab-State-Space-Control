% #########################################################################
% TUHH :: Institute for Control Systems :: Control Lab
% #########################################################################
% Experiment CSTD 1: Observer-based Linear State-feedback Control of TCLab
%
% Copyright Timm Faulwasser and Hamburg University of Technology
% #########################################################################
%
%

%% Task 1.1 Observer-based state feedback design

A = [0.9954 -0.0017 -0.0043 -0.0001;
    0.003 0.9939 0.0026 -0.0129;
    0.0052 -0.0121 0.9705 0.0131;
    0.0022 0.0135 0.009 0.9629];
 
B = [0.0022 0.0013;
    -0.0039 0.0015;
    -0.0041 -0.003;
    0.0029 -0.0026];
 
C = [0.987 -0.1651 0.486 -0.2828;
    0.896 0.4392 0.1921 0.3271];
 
D = [-0.0004 0.0001;
    -0.0025 -0.002];

dt = 1;

% determine controlled equilibrium
T_setpoint = 50;

T_amb      = 23.45;

bar_y = [T_setpoint - T_amb; T_setpoint - T_amb];

% Solve 6x6 system:
% (A-I)*bar_x + B*bar_u = 0
% C*bar_x + D*bar_u = bar_y
I = eye(4);
M = [(A - I), B; C, D];
RS = [0; 0; 0; 0; bar_y];
res = M \ RS;

bar_x = res(1:4);
bar_u = res(5:6);

disp('--- Task 1.1 ---')
disp('bar_x ='); disp(bar_x)
disp('bar_u ='); disp(bar_u)

%% Task 1.2 Verify Controlability
% use ctrb command to get contalability matrix

Ctr1 = B;
Ctr2 = A * Ctr1;
Ctr3 = A * Ctr2;
Ctr4 = A * Ctr3;
controlability_matrix = [Ctr1 Ctr2 Ctr3 Ctr4];

% check the rank and report if it is controlable
control_rank = rank(controlability_matrix);
n = size(A, 1);

if control_rank == n
    disp("rank = 4, controllable")
elseif control_rank < n
    disp("rank less than " + n + ", NOT controllable, number of uncontrollable states = " + string(n - control_rank))
else
    disp("rank > 4, Error")
end

%% Task 1.3 Controller Design
% check eigen values of the plant using eig(A) and use place command for
e = eig(A);
disp('--- Task 1.3 ---')
disp('Open-loop eigenvalues of A:'); disp(e)

% pole placement
%pc = [0.96; 0.94; 0.92; 0.90];
pc = e.^3;

K = place(A, B, pc);
CL_Poles = eig(A - B*K);

disp('Controller gain K ='); disp(K)
disp('Closed-loop poles A-BK ='); disp(CL_Poles)

%% Task 1.4 Verify Observability
% use ctrb command to get observability matrix

Obs1 = C;
Obs2 = Obs1 * A;
Obs3 = Obs2 * A;
Obs4 = Obs3 * A;
observability_matrix = [Obs1; Obs2; Obs3; Obs4];

% check the rank and report if it is observable
obs_Rank = rank(observability_matrix);

if obs_Rank == n
    disp("rank = 4, Obervable")
elseif obs_Rank < n
    disp("rank less than " + n + ", NOT Observable, number of unobservable states = " + string(n - obs_Rank))
else
    disp("rank > 4, Error")
end

%% Task 1.5 Observer Design

%po = [0.88; 0.85; 0.82; 0.80];
po = pc.^2;

L = place(A', C', po)';

disp('--- Task 1.5 ---')
disp('Observer gain L ='); disp(L)
disp('Observer error dynamics eigenvalues (A-LC):')
disp(eig(A - L*C))

%% Task 1.6 Performance evaluation of the control design in simulation

t_control = 300;

u_ctrl_sim_data = zeros(2, t_control);
y_ctrl_sim_data = zeros(2, t_control);
observer_error_sim = zeros(2, t_control);
tracking_error_sim = zeros(2, t_control);

n_x = 4;
hat_x = zeros(n_x, 1);

mdl_x = zeros(n_x, 1);

    for counter = 1:t_control

        % calculate u_ctrl and consider saturation
        u_ctrl =  bar_u - K * (hat_x - bar_x);
        u_ctrl = max(0, min(100, u_ctrl));

        % define noise function
        Noise_cov = 0.01 * eye(2);
        y_signal = C * mdl_x + D * u_ctrl + chol(Noise_cov) * randn(2, 1); % consider measurement noise

        % Record output-based observer and tracking errors for simulation.
        % These quantities correspond to the error plots reported in the experiment.
        hat_y = C * hat_x + D * u_ctrl;
        observer_error_sim(:, counter) = y_signal - hat_y;
        tracking_error_sim(:, counter) = bar_y - y_signal;
    

        % apply u_ctrl to observer
        hat_x = A * hat_x + B * u_ctrl + L * (y_signal - hat_y);
        
        % apply u_ctrl to model
        mdl_x = A * mdl_x + B * u_ctrl;
    
        u_ctrl_sim_data(:, counter) = u_ctrl;
        y_ctrl_sim_data(:, counter) = y_signal;
    
    end

%% plot simulation

t = 1:t_control;

figure(100);
clf;

sgtitle('Simulation Results: Output Temperatures and Control Inputs');

subplot(2,2,1);
plot(t, y_ctrl_sim_data(1,:) + T_amb, 'LineWidth', 1.5);
hold on;
yline(T_setpoint, '--r', 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Temperature [^\circC]');
title('Heater 1 Output Temperature');
legend('T_1', 'Setpoint', 'Location', 'best');

subplot(2,2,2);
plot(t, u_ctrl_sim_data(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Input [%]');
title('Heater 1 Control Input');
legend('u_1', 'Location', 'best');

subplot(2,2,3);
plot(t, y_ctrl_sim_data(2,:) + T_amb, 'LineWidth', 1.5);
hold on;
yline(T_setpoint, '--r', 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Temperature [^\circC]');
title('Heater 2 Output Temperature');
legend('T_2', 'Setpoint', 'Location', 'best');

subplot(2,2,4);
plot(t, u_ctrl_sim_data(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Input [%]');
title('Heater 2 Control Input');
legend('u_2', 'Location', 'best');


figure(101);
clf;

sgtitle('Simulation Results: Observer Error and Tracking Error');

subplot(2,2,1);
plot(t, observer_error_sim(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Observer Error - Output 1');
legend('e_{obs,1}', 'Location', 'best');

subplot(2,2,2);
plot(t, tracking_error_sim(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Tracking Error - Output 1');
legend('e_{trk,1}', 'Location', 'best');

subplot(2,2,3);
plot(t, observer_error_sim(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Observer Error - Output 2');
legend('e_{obs,2}', 'Location', 'best');

subplot(2,2,4);
plot(t, tracking_error_sim(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Tracking Error - Output 2');
legend('e_{trk,2}', 'Location', 'best');

%% Task 1.7 Performance evaluation of the control design on the real system

u_ctrl_data = zeros(2, t_control);
y_ctrl_data = zeros(2, t_control);
observer_error = zeros(2, t_control);
tracking_error = zeros(2, t_control);

hat_x = zeros(n_x,1);

% Connect to TCLab
lab = tclab;

    for counter = 1:t_control    

        tStart = tic;

        % measure the current output
        y_abs = [lab.T1; lab.T2];
        y_signal = y_abs - T_amb;
        %* ones(2,1) - (C * bar_x + D * bar_u);

        % calculate u_ctrl considering saturation
        u_ctrl = bar_u - K * (hat_x - bar_x);
        u_ctrl = max(0, min(100, u_ctrl));
    
        % turn off hearters and report warning message if the temperatrue of
        % any sensor exceeds 70 C
        if y_abs(1) > 70 || y_abs(2) > 70
            u_ctrl = [0; 0];
            lab.Q1(0);
            lab.Q2(0);
            warning('Temerature exceeded 70 C! Heaters turned off at step %d.', counter);
        end
      
        % apply u_ctrl in observer dynamics
        hat_y = C * hat_x + D * u_ctrl;
        hat_x = A * hat_x + B * u_ctrl + L * (y_signal - hat_y);
   
        % apply u_ctrl to TClab
        lab.Q1(u_ctrl(1));
        lab.Q2(u_ctrl(2));

        u_ctrl_data(:, counter) = u_ctrl;
        y_ctrl_data(:, counter) = y_signal;
        observer_error(:, counter) = y_signal - hat_y;

        tracking_error(:, counter) = bar_y - y_signal; %deviation from setpoint
    
    
        % let the LED blink for 0.1 s
        lab.blink(0.1);

        disp(['Counter ', int2str(counter)]);
        disp(['U= ', num2str(u_ctrl.')]);
        disp(['T= ', num2str(y_abs.')]);

        while (toc(tStart) < dt)
        
        end
    
    end

lab.Q1(0);
lab.Q2(0);
disp(lab.off());
clear lab;

%% plot real-time control

t = 1:t_control;

figure(200);
clf;

sgtitle('Real-System Results: Output Temperatures and Control Inputs');

subplot(2,2,1);
plot(t, y_ctrl_data(1,:) + T_amb, 'LineWidth', 1.5);
hold on;
yline(T_setpoint, '--r', 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Temperature [^\circC]');
title('Heater 1 Output Temperature');
legend('T_1', 'Setpoint', 'Location', 'best');

subplot(2,2,2);
plot(t, u_ctrl_data(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Input [%]');
title('Heater 1 Control Input');
legend('u_1', 'Location', 'best');

subplot(2,2,3);
plot(t, y_ctrl_data(2,:) + T_amb, 'LineWidth', 1.5);
hold on;
yline(T_setpoint, '--r', 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Temperature [^\circC]');
title('Heater 2 Output Temperature');
legend('T_2', 'Setpoint', 'Location', 'best');

subplot(2,2,4);
plot(t, u_ctrl_data(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Input [%]');
title('Heater 2 Control Input');
legend('u_2', 'Location', 'best');

figure(201);
clf;

sgtitle('Real-System Results: Observer Error and Tracking Error');

subplot(2,2,1);
plot(t, observer_error(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Observer Error - Output 1');
legend('e_{obs,1}', 'Location', 'best');

subplot(2,2,2);
plot(t, tracking_error(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Tracking Error - Output 1');
legend('e_{trk,1}', 'Location', 'best');

subplot(2,2,3);
plot(t, observer_error(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Observer Error - Output 2');
legend('e_{obs,2}', 'Location', 'best');

subplot(2,2,4);
plot(t, tracking_error(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time step k');
ylabel('Error');
title('Tracking Error - Output 2');
legend('e_{trk,2}', 'Location', 'best');