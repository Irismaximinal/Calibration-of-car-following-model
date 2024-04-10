function [f] = fitness(x,dt,dv,dx,v,x_l,v_l)
% x: parameter vector
% x(1): desired time_headway T
% x(2): desired speed: v0
% x(3): acceleration exponential coeffcient:  δ
% x(4): maximun acceleration: α
% x(5): comfortable deceleration: 	β
% x(6): minimum gap: s0
% formulation of IDM 
% dv/dt = α*( 1-(v/v0)^δ-( ( s0+vT+vΔv/sqrt(|4αβ|) )/ (dx-l) )^2 )

gama = 0.5;
len = length(dv);
car_len = 6;

% Simulation state vector
dx_ = dx(1)*ones(size(dx));
dv_ = dv(1)*ones(size(dv));
v_ = v(1)*ones(size(v));
x_ = (x_l(1)-dx(1))*ones(size(x_l));

% Trajectory simulation
for k = 1:len-1
    
    % Acceleration calculation
    e_dis = x(6)+v_(k)*x(1) - v_(k)*dv_(k)/sqrt(-4*x(4)*x(5));
    a_k = x(4)*( 1- ( v_(k)/x(2) )^x(3) - ( e_dis/(dx_(k)-car_len) )^2);

     % Position and velocity updating
    v_(k+1) = v_(k)+a_k*dt;
    x_(k+1) = x_(k)+v_(k)*dt+0.5*a_k*dt^2;
    dv_(k+1) = v_l(k+1)-v_(k+1);
    dx_(k+1) = x_l(k+1)-x_(k+1);
    
end

% Finess function: to mininize the speed and gap error
% Including:
% percentage error of related speed: (dv-dv_)./dv
% percentage error of gap: (dx-dx_)./dx
% Root mean square percentage error (RMSPE)

RMSPE_dv = mean(((dv-dv_)./(0.01+dv)).^2);
RMSPE_dx = mean(((dx-dx_)./(0.01+dx)).^2);

f = gama*RMSPE_dx + (1-gama)*RMSPE_dv;
