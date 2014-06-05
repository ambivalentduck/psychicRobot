function dq=armdynamicsShadMuss(t,q)

global measuredVals measuredTime kpgain

lq=length(q);

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
theta_desired=interped(1:2)';
omega_desired=interped(3:4)';
alpha_desired=interped(5:6)';
torque_outside=interped(7:8)';

% Add feedback forces
theta_real=q(1:lq/2);
omega_real=q(lq/2+1:end);

% Compute alpha to torque relationship
[D_real,C_real]=computeDC(theta_real,omega_real);
[D_expected,C_expected]=computeDC(theta_desired,omega_desired);

kp=[15 6;6 16]*kpgain;
kd=[2.3 .09; .09 2.4]*kpgain;

torque_fb=kp*(theta_real-theta_desired) + kd*(omega_real-omega_desired);

%Update the change in desired state
dq=[omega_real;
    D_real\(D_expected*alpha_desired+C_expected-C_real-torque_fb-torque_outside);];  %If torque_fb and torque_outside=0, and c_real ~ c_expected, alpha = alpha desired.