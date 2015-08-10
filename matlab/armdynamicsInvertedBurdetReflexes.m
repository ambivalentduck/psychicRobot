function [dqi,torque_fb,torque_inertial]=armdynamicsInvertedBurdetReflexes(t,qi)

global measuredVals measuredTime errorVals errorTime kpgain massgain reflexcontrib

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
theta_real=interped(1:2)';
omega_real=interped(3:4)';
alpha_real=interped(5:6)';
torque_outside=interped(7:8)';

interped2=twoNearestNeighbor(errorVals,errorTime,t-.06);
reflexE=interped2(1:2)';
reflexV=interped2(3:4)';

% Add feedback forces
theta_desired=qi(1:2);
omega_desired=qi(3:4);

% Compute alpha to torque relationship
[D_real,C_real]=computeDC(theta_real,omega_real);
[D_expected,C_expected]=computeDC(theta_desired,omega_desired);

kp0=[10.8 2.83; 2.51 8.67];
joint_torques=abs(D_real*alpha_real+C_real+torque_outside);
kp1=[3.18 2.15; 2.34 6.18];
kp=kp0+kp1*diag(joint_torques);

kp=kpgain*kp;

%should be (kp/#50#)*(reflexE+2*reflexV)
torque_fb=(1-reflexcontrib)*kp*((theta_real-theta_desired) + (1/12)*(omega_real-omega_desired))+(reflexcontrib*kp)*(reflexE+2*reflexV);


% Update the change in desired state
dqi=[omega_desired;
    D_expected\(D_real*alpha_real+(torque_fb+torque_outside)/massgain+C_real-C_expected);];  %If torque_fb and torque_outside=0, and c_real ~ c_expected, alpha = alpha desired.

torque_inertial=D_expected*dqi(3:4)+C_expected;
end
