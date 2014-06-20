function [dqi,torque_outside]=armdynamics_inverted(t,qi)

global K1diag K2diag Koffdiag measuredVals measuredTime

lqi=length(qi);

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
theta_real=interped(1:2)'; 
omega_real=interped(3:4)';
alpha_real=interped(5:6)';
torque_outside=interped(7:8)';

%Add feedback forces
theta_desired=qi(1:2);
omega_desired=qi(3:4);

%Compute alpha to torque relationship, eq. 7.87 in Spong's Robot Control and Modeling: pg 262
[D_real,C_real]=computeDC(theta_real,omega_real);
[D_expected,C_expected]=computeDC(theta_desired,omega_desired);

kp=[K1diag Koffdiag; Koffdiag K2diag];
kd=[2.3 .09; .09 2.4];

torque_fb=kd*(omega_desired-omega_real)+kp*(theta_desired-theta_real);

%Update the change in desired state
dqi=[omega_desired;
    D_expected\(D_real*alpha_real-torque_fb+torque_outside+C_real-C_expected);];  %If torque_fb and torque_outside=0, and c_real ~ c_expected, alpha = alpha desired.
end
