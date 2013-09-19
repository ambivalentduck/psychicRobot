function [dq,inertial]=armdynamicsBurdetReflexes(t,q)

global measuredVals measuredTime errorVals errorTime inertialTorque kpgain

lq=length(q);

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
theta_desired=interped(1:2)';
omega_desired=interped(3:4)';
alpha_desired=interped(5:6)';
torque_outside=interped(7:8)';

interped2=twoNearestNeighbor(errorVals,errorTime,t-.06);
reflexE=interped2(1:2)';
reflexV=interped2(3:4)';

% Add feedback forces
theta_real=q(1:lq/2);
omega_real=q(lq/2+1:end);

% Compute alpha to torque relationship
[D_real,C_real]=computeDC(theta_real,omega_real);
[D_expected,C_expected]=computeDC(theta_desired,omega_desired);

kp0=[10.8 2.83; 2.51 8.67];
% In inverse version, the inertial and outside torques have already been measured
% joint_torques=abs(D_real*alpha_real+C_real+torque_outside);
% Instead introduce a tiny lag in the stiffness by relying on .
joint_torques=abs(inertialTorque+torque_outside);
kp1=[3.18*joint_torques(1) 2.15*joint_torques(2); 2.34*joint_torques(1) 6.18*joint_torques(2)];
kp=kpgain*(kp0+kp1);

torque_fb=kp*((theta_real-theta_desired) + (1/12)*(omega_real-omega_desired))+(kp/50)*(reflexE+2*reflexV);

%Update the change in desired state
dq=[omega_real;
    D_real\(D_expected*alpha_desired+C_expected-C_real-torque_fb-torque_outside);];  %If torque_fb and torque_outside=0, and c_real ~ c_expected, alpha = alpha desired.

if nargout>1
    inertial=D_real*dq(3:4)+C_real;
end
