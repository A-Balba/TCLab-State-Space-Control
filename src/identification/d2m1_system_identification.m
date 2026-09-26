% #########################################################################
% TUHH :: Institute for Control Systems :: Control Lab
% #########################################################################
% Experiment D2M 1: Linear and Nonlinear Model Identification of TCLab
%
% Copyright Timm Faulwasser and Hamburg University of Technology
% #########################################################################
%
%
% 'help command' into the MATLAB window.

%% Task 1.1: Data Loading and Investigation
% Purpose: Load provided identification/validation datasets and visualize
%          the signals to become familiar with sampling, scales, and trends.
clc
clear

% Load experiment data (provided in data.mat).
% Expected variables include: dt_data, n_samples, n_val,
%                             y_data (2 x n_samples), u_data (2 x n_samples),
%                             y_val  (2 x n_val),     u_val  (2 x n_val),
%                             T_amb_data (scalar ambient temperature in K).
load("data.mat")
%input-output data curve (7 hours !)
%validation data is 1 h 20 min. in comparison !
%Data Analysis.
% variables in training data (states)
% T1: plant 1 temperature
% Q1: Heater 1 power percentage
% t2: plant 2 temperature
% Q2: Heater 2 power percentage
% T10 = 299.885 k =  @ Time=0
% T1 ambient= 26.7 C
% T20 = T10 = 299.4643 k
% T2 ambient = 26.3143 C
%Sampling time = 20 s


% Plot identification dataset (about 7 hours)
figure(100)
sgtitle("7 Hours of Measured Training Data")
subplot(2,2,1)
plot(1:dt_data:n_samples*dt_data, y_data(1,:))
grid("on")
title('T_1 (K)')

subplot(2,2,2)
plot(1:dt_data:n_samples*dt_data, y_data(2,:))
grid("on")
title('T_2 (K)')

subplot(2,2,3)
plot(1:dt_data:n_samples*dt_data, u_data(1,:))
grid("on")
title('Q_1 (%)')
xlabel('Time (s)')

subplot(2,2,4)
plot(1:dt_data:n_samples*dt_data, u_data(2,:))
grid("on")
title('Q_2 (%)')
xlabel('Time (s)')

% Plot validation dataset (about 40 minutes)
figure(101)
sgtitle("80 Minutes of Measured Validation Data")
subplot(2,2,1)
plot(1:dt_data:n_val*dt_data, y_val(1,:))
grid("on")
title('T_1 (K)')

subplot(2,2,2)
plot(1:dt_data:n_val*dt_data, y_val(2,:))
grid("on")
title('T_2 (K)')

subplot(2,2,3)
plot(1:dt_data:n_val*dt_data, u_val(1,:))
grid("on")
title('Q_1 (%)')
xlabel('Time (s)')

subplot(2,2,4)
plot(1:dt_data:n_val*dt_data, u_val(2,:))
grid("on")
title('Q_2 (%)')
xlabel('Time (s)')

%% Task 1.2: Ambient Temperature Calibration
% Purpose: Re-estimate the ambient temperature and noise covariance directly
%          from the hardware at the time of the experiment. This helps align
%          the dataset with current lab conditions.
clc
lab = tclab ();

n  = 100;     % number of ambient measurements
dt = 1;       % seconds between measurements

T_meas_array = zeros(2, n);

for counter = 1:n
    % Simple timing guard to achieve ~1 s sampling
    tStart = tic;
    % Read current temperatures (in degC). lab.T1 and lab.T2 are in degC.
    T_meas_array(:, counter) = [lab.T1; lab.T2];
    disp(['T1 = ', num2str(T_meas_array(1, counter)), '  T2 = ', num2str(T_meas_array(2, counter))]);
    while (toc(tStart) < dt)
        % busy-wait until 1 s elapses
    end
end

% Turn off heaters and close connection
disp(lab.off());
%% 
% Convert to Kelvin and compute ambient estimate
%initialize T1_avg and T2_avg
T1_avg=0;
T2_avg=0;
    for i=1:n
    T1_avg = T1_avg + T_meas_array (1,i);
    T2_avg = T2_avg + T_meas_array (2,i);
    end
T1_avg = T1_avg/(n);
T2_avg = T2_avg/(n);

T_amb_experiment = (T1_avg + T2_avg)/2 ;
T_amb_experiment = T_amb_experiment+273.15;

% Empirical covariance of sensor noise (degC^2). Keep for reference.
Noise_cov = cov(T_meas_array.');

disp(['Experiment ambient temperature is ', num2str(T_amb_experiment), ' K (', num2str(T_amb_experiment - 273.15), ' °C).'])
clear dt counter lab tStart T_meas_array n
%Tamb @ experiment = 301.8282
%% Task 1.3: Data Calibration
% Purpose: Align the provided data with the newly measured ambient
%          temperature so that absolute levels are consistent.
T_amb_diff = T_amb_experiment-T_amb_data;
%   make two arrays for calibration of data and validation
T_amb_diff_array_data = T_amb_diff * ones(size(y_data));
T_amb_diff_array_val  = T_amb_diff * ones(size(y_val));
% calibration
y_data = y_data + T_amb_diff_array_data;
y_val  = y_val + T_amb_diff_array_val;

%% Task 2.1: Linear Subspace Model Identification

L =200;
%   The embedding dimension (block size), which should be greater 
%   than the expected system order. 
[sigma, ssfun] = moesp(y_data', u_data', L);

bar((sigma'))

set(gca,'yscal','log');

%% Task 2.2: Performance Evaluation of Identified Models

% Propose the smallest system order here
% 1 is chosen
n_x_lb = 1;
% Highest chosen system order
% 6 is chosen
n_x_ub = n_x_lb + 5;

y_mse = zeros(1, n_x_ub - n_x_lb + 1);


for n_x_choice = n_x_lb:n_x_ub

    [An,Bn,Cn,Dn] = ssfun(n_x_choice);
    O_L = observability_matrix (An, Cn,10);
    H = toeplitz_tril (An, Bn, Cn, Dn, 10);
    % Get the initial condition with respect to the estimated system matrices
    x_0= pinv(O_L)*(y_val(1:20)'-H*(u_val(1:20)'));

    % Get simulated output trajectories given the validation input
    % u_val_data, the initial condition x_0, and the system matrices An, Bn,
    % Cn, and Dn
    
    x_sim(:,1) = x_0;

    for i = 2:n_val
        x_sim(:,i) = An*x_sim(:,i-1)+Bn*u_val(:,i-1);
    end
    y_sim = Cn*x_sim+Dn*u_val;
    clear x_sim

    % Calculate the mean square error of your simulated trajectory compared
    % to the validation trajectory
    y_mse(n_x_choice - n_x_lb + 1) = mean(mean((y_sim'-y_val').^2));%immse(y_sim,y_val);
end

% Plot the mean square errors
figure
plot([n_x_lb:n_x_ub], y_mse)
xlabel('System Order');
ylabel('MSE')
title('Mean Square Error');
grid on
y_mse

%% Task 2.3: Choosing the Best Linear Model Order and Evaluating Its Performance
% Selected model order used for validation.
n_x = 3;
[An,Bn,Cn,Dn] = ssfun(n_x);
O_L = observability_matrix (An, Cn,10);
H = toeplitz_tril (An, Bn, Cn, Dn, 10);
x_0= pinv(O_L)*(y_val(1:20)'-H*(u_val(1:20)'));

clear x_sim
x_sim(:,1) = x_0;

for i = 2:n_val

    x_sim(:,i) = An*x_sim(:,i-1) + Bn*u_val(:,i-1);
end

y_est_linear_val =Cn*x_sim+Dn*u_val;
y_final_error = immse(y_est_linear_val,y_val)
% Simple comparison plots
figure(200)

subplot(2,2,1)
plot(1:dt_data:(n_val-1)*dt_data, (y_val(1,2:n_val)' - y_val(1,1:n_val-1)') / dt_data)
hold all
plot(1:dt_data:(n_val-1)*dt_data, (y_est_linear_val(1,2:n_val)' - y_est_linear_val(1,1:n_val-1)') / dt_data)
title('dT_1/dt (approx.)')
legend('Validation','Estimate')

subplot(2,2,2)
plot(1:dt_data:(n_val)*dt_data, y_val(1,1:end))
hold all
plot(1:dt_data:(n_val)*dt_data, y_est_linear_val(1,1:end))
title('T_1 (K)')
legend('Validation','Estimate')

subplot(2,2,3)
plot(1:dt_data:(n_val-1)*dt_data, (y_val(2,2:n_val)' - y_val(2,1:n_val-1)') / dt_data)
hold all
plot(1:dt_data:(n_val-1)*dt_data, (y_est_linear_val(2,2:n_val)' - y_est_linear_val(2,1:n_val-1)') / dt_data)
xlabel('Time (s)')
title('dT_2/dt (approx.)')
legend('Validation','Estimate')

subplot(2,2,4)
plot(1:dt_data:(n_val)*dt_data, y_val(2,1:end))
hold all
plot(1:dt_data:(n_val)*dt_data, y_est_linear_val(2,1:end))
xlabel('Time (s)')
title('T_2 (K)')
legend('Validation','Estimate')

%% Task 3.1: Data Organization for Kernel Regression
% Purpose: Construct the regressor vector x_k that stacks m past outputs and
%          the current inputs, i.e.,
%          x_k = [y(k-1); y(k-2); ...; y(k-m); u(k-1)], with m = 2.

%T_amb: assume 300 kelvin, in Lab it will be measured
%T_amb = T_amb_data+T_amb_diff;      % Store ambient temperature for later use in kernel
T_amb = T_amb_experiment
m = 2;                              % Number of output lags, for example 2

%Kvecotr has last two data points since m=2 (subject to change)
%T1 : 2 rows 
%T2 : 2 rows
%Q1 : 1 coloumn
%Q2 : 1 coloumn
%number of coloumns is the sample size minus m
x = zeros (2*(m+1),n_samples-m);
for k=m+1:n_samples
x(:,k-m) = [y_data(:,k-1);
            y_data(:,k-2);
            u_data(:,k-1)];
end
%% Task 3.2: Kernel Function
% Implemented in kernel.m as a degree-2 polynomial kernel with ambient centering.

%% Task 3.3: Gram Matrix Function
% Implemented in gram_matrix.m using the kernel defined above.

%% Task 3.4: Kernel Regression Model Identification
% Purpose: Compute kernel weights (alpha) that map regressors to next-step
%          outputs via regularized least squares in the RKHS.
rho = 0.0001;                                  % regularization
v   = y_data(:, m+1:end).';                                  % targets: y(k)
G   = gram_matrix(x, T_amb);                % N x N
alpha = (G + rho*eye(size(G))) \ v;                                % N x 2 (two outputs)

%% Task 3.5: Prediction Function
% Implemented in f_next.m.

%% Task 3.6: Validation Using Validation Data
% Purpose: Perform one-step-ahead simulation on the validation dataset.

x_val = zeros (2*(m+1), n_val); % NOTE: uses same sizing as original
y_est_nonlinear_val = zeros (2, n_val);
y_est_nonlinear_val(:,1:m) = y_val(:,1:m); % seed with first m measurements

for counter = m+1:n_val

        % Shift window and append new measurement and known input
    x_val(:, counter) = [y_val(:,counter-1);
                        y_val(:,counter-2);
                        u_val(:,counter-1)];
    
    % One-step prediction using identified model
    y_est_nonlinear_val(:, counter) = f_next(x_val(:,counter), x, T_amb, alpha);

end

% Simple comparison plots
figure(300)

subplot(2,2,1)
plot(1:dt_data:(n_val-1)*dt_data, (y_val(1,2:n_val)' - y_val(1,1:n_val-1)') / dt_data)
hold all
plot(1:dt_data:(n_val-1)*dt_data, (y_est_nonlinear_val(1,2:n_val)' - y_est_nonlinear_val(1,1:n_val-1)') / dt_data)
title('dT_1/dt (approx.)')
legend('Validation','Estimate')

subplot(2,2,2)
plot(1:dt_data:(n_val)*dt_data, y_val(1,1:end))
hold all
plot(1:dt_data:(n_val)*dt_data, y_est_nonlinear_val(1,1:end))
title('T_1 (K)')
legend('Validation','Estimate')

subplot(2,2,3)
plot(1:dt_data:(n_val-1)*dt_data, (y_val(2,2:n_val)' - y_val(2,1:n_val-1)') / dt_data)
hold all
plot(1:dt_data:(n_val-1)*dt_data, (y_est_nonlinear_val(2,2:n_val)' - y_est_nonlinear_val(2,1:n_val-1)') / dt_data)
xlabel('Time (s)')
title('dT_2/dt (approx.)')
legend('Validation','Estimate')

subplot(2,2,4)
plot(1:dt_data:(n_val)*dt_data, y_val(2,1:end))
hold all
plot(1:dt_data:(n_val)*dt_data, y_est_nonlinear_val(2,1:end))
xlabel('Time (s)')
title('T_2 (K)')
legend('Validation','Estimate')

%% Task 4.1: MSE Comparison of Linear and Nonlinear Models (Validation Data)

% --- Linear model (already computed as y_est_linear_val) ---
y_mse_linear_val    = mean(mean((y_est_linear_val'    - y_val').^2));

% --- Nonlinear model (kernel-based, from Task 3.6) ---
y_mse_nonlinear_val = mean(mean((y_est_nonlinear_val' - y_val').^2));

% Display results
fprintf('\n=============================================\n');
fprintf('   VALIDATION DATA PERFORMANCE (MSE)\n');
fprintf('=============================================\n');

fprintf('Nonlinear Model MSE (Validation): %.6f\n', y_mse_nonlinear_val);
fprintf('Linear Model MSE    (Validation): %.6f\n', y_mse_linear_val);

% Relative improvement
improvement_val = (y_mse_linear_val - y_mse_nonlinear_val) / y_mse_linear_val * 100;

fprintf('---------------------------------------------\n');
fprintf('Relative Improvement (Nonlinear vs Linear): %.2f %%\n', improvement_val);

% Conclusion
if y_mse_nonlinear_val < y_mse_linear_val
    fprintf('Conclusion: Nonlinear model performs better on validation data.\n');
elseif y_mse_nonlinear_val > y_mse_linear_val
    fprintf('Conclusion: Linear model performs better on validation data.\n');
else
    fprintf('Conclusion: Both models perform equally on validation data.\n');
end

fprintf('=============================================\n\n');

