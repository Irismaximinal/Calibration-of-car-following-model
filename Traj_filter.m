function [s_,v_] = Traj_filter(t,dt,s,v)
% to obtain the smoothed trajectory by some filters
% Finished filter: Kalman filter
% Step of kalman filter including:
% Step 1: Prediction
% get the new system state x and prior covariance maxtrix P
% x_pred = A*x + B*u
% p_pred = A*P*A' + Q
% Step 2: Updating
% get the new value of Kalman gain K, optimal estimated state x, 
% posteriori covariance maxtrix P
% K = P_pred * H' * inv(H * P_pred * H' + R)
% x = x_pred + K * (measurements(i, :)' - H * x_pred);
% P = (eye(2) - K * H) * P_pred;

%%  Modeling and Parameter setting

% Noisy matrix
Q = 0.1*eye(2);                    % process noise, the matrix Q is a hyperparameter
R = [1, 0; 0, 1];                      % measure noise, the matrix R is a hyperparameter

% Predicted matrix
% Linear state space model/ bicycle model
% Obtained by System identification
A = [1, dt; 0, 1];                  % state transition matrix
B = [0.5*dt^2; dt];             % control input matrix
u = 0;                                   % control input/ set as zero for unknown value 

% Measurement matrix
H = eye(2);                         % measurement matrix

% Initial covariance matrix
P = eye(2);

% Measurement data
x = [s,v];

%%  Kalman estimation
x0 = [s(1); v(1)];                  % initial state
x_ = x0'.*ones(size(x));      % estimated data
x_k = [s(2); v(2)];                 % initial state in kalman filter, x_k = x(2) if x0 is accurate

% loop for state estimation
for k = 2:size(t)
    % Prediction
    x_pred = A*x_k + B*u;
    P_pred = A*P*A' + Q;

    % Updating
    K = P_pred*H'*inv(H*P_pred*H' + R);
    x_k = x_pred + K*(x(k, :)' - H*x_pred);
    P = (eye(2) - K*H)*P_pred;

    % Optimal estimation storage
    x_(k,:) = x_k';
end

s_ = x_(:,1);
v_ = x_(:,2);

%% filter visualization
% figure('color','w');
% [ax,hLine1,hLine2]=plotyy(t, s, t, v);
% hLine1.LineStyle = '--';
% hLine2.LineStyle = '--';
% hold on
% [ax,hLine3,hLine4]=plotyy(t, s_, t, v_);
% hLine3.LineStyle = '-';
% hLine4.LineStyle = '-';
end