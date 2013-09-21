function [dqi]=armdynamicsInvertedBurdetReflexes(t,qi)

global measuredVals measuredTime errorVals errorTime kp0 kp1 kpgain

lqi=length(qi);

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
theta_real=interped(1:2)';
omega_real=interped(3:4)';
alpha_real=interped(5:6)';
torque_outside=interped(7:8)';

interped2=twoNearestNeighbor(errorVals,errorTime,t-.06);
reflexE=interped2(1:2)';
reflexV=interped2(3:4)';

% Add feedback forces
theta_desired=qi(1:lqi/2);
omega_desired=qi(lqi/2+1:end);

% Compute alpha to torque relationship
[D_real,C_real]=computeDC(theta_real,omega_real);
[D_expected,C_expected]=computeDC(theta_desired,omega_desired);


joint_torques=abs(D_real*alpha_real+C_real+torque_outside);
kpT=kp1*[joint_torques(1) 0; 0 joint_torques(2)];
kp=kpgain*(kp0+kpT);

torque_fb=kp*((theta_real-theta_desired) + (1/12)*(omega_real-omega_desired))+(kp/50)*(reflexE+2*reflexV);

% Update the change in desired state
dqi=[omega_desired;
    D_expected\(D_real*alpha_real+torque_fb+torque_outside+C_real-C_expected);];  %If torque_fb and torque_outside=0, and c_real ~ c_expected, alpha = alpha desired.
end
