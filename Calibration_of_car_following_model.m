%---------Car-following model calibration----------%
%--------- Author: Changxin Wan---------%
clc;clear;
dbstop if error;
%% Data preparation
%----------data loading----------%
% Loading the raw data from you file
% Including: time, position, velocity of leader and follower vehicle

path = 'C:\Users\irism\Desktop\Car_following_data\';        % position where you data storage
filename = 'car_following_dataset.csv';                                    % name and type of the data

dataset = importdata([path,filename]);

time = dataset.data(:,1);   
% position_lead = dataset.data(:,2);                        
% velocity_lead = dataset.data(:,3);                        
% position_follower = dataset.data(:,4);                 
% velocity_follower = dataset.data(:,5); 

dt = time(2)-time(1);                                              % time step 

position_lead = dataset.data(:,2)+normrnd(0,0.5,size(time));                        
velocity_lead = dataset.data(:,3)+normrnd(0,0.1,size(time));                        
position_follower = dataset.data(:,4)+normrnd(0,0.5,size(time));                 
velocity_follower = dataset.data(:,5)+normrnd(0,0.05,size(time));                  

%----------data filter/ Unfinished----------%
% The raw data may have some noise 
% Using the filter to smooth the trajectory
% Filter:  Kalman filter

[x_l,v_l] = Traj_filter(time,dt, position_lead, velocity_lead);
[x_f,v_f] = Traj_filter(time,dt, position_follower, velocity_follower);

%----------car_following data transformation----------%
% Transfer the smoothed trajectroy to the data used in car following model
% Including: time, position, velocity of leader and follower vehicle
dx = x_l-x_f;                               % bump to bump distance
dv = v_l-v_f;                               % related speed
v=v_f;                                          % current speed

len = length(dx);                       % Total number of time step
len_cal = round(len*0.8);

%% Model calibration
%----------car following parameter----------%
% parameter needed to be calibrated
% car-following model: IDM
% equation of IDM: 
% dv/dt = alpha*( 1- (v/v0)^sita - ( (s0+v*h_safe+v*dv/sqrt(alpha*beta))/(dx-l) )^2 )
% Parameter including:
% desired time_headway : h_desire, unit :s
% desired speed: v0, unit : m/s
% acceleration exponential coeffcient: sita
% maximun acceleration: alpha, unit : m/s^2
% comfortable deceleration: beta, unit : m/s^2
% minimum gap: s0, unit : m


% lower and upper bound of parameters
h_safe = [1.1, 2.6];
v0 = [15.0, 33.3];
sita = [-10.0, 10.0];
alpha = [0.5,3.0];
beta = [-2.0,-0.5];
s0 = [1.0, 3.0];

X = [h_safe;v0;sita;alpha;beta;s0];
NumVars = length(X(:,1));                 % number of calibration parameters
X_L = X(:,1);                                           % lower bound of parameters
X_U = X(:,2);                                          % upper bound of parameters

%----------Genetic algorithm parameter----------%
% hyperparameter for ga
% Selection method: roulette
% Elite number: 3
% Crossover method: twopoint
% Crossover possibility: 0.8
% Mutation possibility: 0.05

options = optimoptions('ga',...
    'SelectionFcn',{@selectionroulette},...
    'EliteCount', 3,...
    'CrossoverFcn','crossovertwopoint',...
    'CrossoverFraction',0.8,...
    'MutationFcn', {@mutationadaptfeasible,0.05},...
    'PlotFcn', @gaplotbestfun);

dx_cal = dx(1:len_cal);                 % bump to bump distance for calibration dataset
dv_cal = dv(1:len_cal);                 % related speed for calibration dataset
v_cal = v(1:len_cal);                      % current speed for calibration dataset
x_cal = x_f(1:len_cal);                   % current position for calibration dataset
x_l_cal = x_l(1:len_cal);                 % leader position for calibration dataset
v_l_cal = v_l(1:len_cal);                 % leader speed for calibration dataset

[X,fval,exitflag,output,population,scores] =...
    ga(@(x) fitness(x,dt,dv_cal,dx_cal,v_cal,x_l_cal,v_l_cal), ...
            NumVars,[],[],[],[],X_L,X_U,[],[],options);

disp('The optimal value of parameters');
fprintf(['h_safe:%f\n' ...
            'v0: %f\n' ...
            'sita: %f\n' ...
            'alpha: %f\n' ...
            'beta: %f\n' ...
            's0: %f\n'],...
            X(1), X(2), X(3), X(4), X(5), X(6));

%% Model validation
% Using the calibrated CF model to simulate the remain trajectroy


car_len = 6;

dx_val = dx(len_cal:end);
dv_val = dv(len_cal:end);
v_val = v(len_cal:end);
x_val = x_f(len_cal:end);
x_l_val = x_l(len_cal:end);
v_l_val = v_l(len_cal:end);

% Simulation state vector
dx_val_ = dx_val(1)*ones(size(dx_val));
dv_val_ = dv_val(1)*ones(size(dv_val));
v_val_ = v_val(1)*ones(size(v_val));
x_val_ = (x_l_val(1)-dx_val(1))*ones(size(x_l_val));

% Trajectory simulation
for k = 1:(len-len_cal)
    
    % Acceleration calculation
    e_dis = X(6)+v_val_(k)*X(1) - v_val_(k)*dv_val_(k)/sqrt(-X(4)*X(5));
    a_k = X(4)*( 1- ( v_val_(k)/X(2) )^X(3) - ( e_dis/(dx_val_(k)-car_len) )^2);

    % Position and velocity updating
    % Using bicycle model:
    % v_(k+1) = v_(k) + a_k*dt
    % x_(k+1) = x_k + v_(k)*dt + 0.5*a_k*dt^2

    v_val_(k+1) = v_val_(k)+a_k*dt;
    x_val_(k+1) = x_val_(k)+v_val_(k)*dt+0.5*a_k*dt^2;
    dv_val_(k+1) = v_l_val(k+1)-v_val_(k+1);
    dx_val_(k+1) = x_l_val(k+1)-x_val_(k+1);
    
end

figure('color','w');
hold on;
plot(time(1:len_cal),x_cal,'-b');
plot(time(len_cal:end),x_val,'--b');
plot(time(len_cal:end),x_val_,'--r');
legend('Previous trajectory','Raw data','IDM');


figure('color','w')
hold on;
plot(time(1:len_cal),v_cal,'-b');
plot(time(len_cal:end),v_val,'--b');
plot(time(len_cal:end),v_val_,'--r');
legend('Previous trajectory','Raw data','IDM');
