function [D,C]=computeDC(theta, omega)

global l1 lc1 lc2 m1 m2 I1 I2

%Compute alpha to torque relationship, eq. 7.87 in Spong's Robot Control
%and Modeling: pg 262

c2=cos(theta(2));
d11=I1+I2+m1*lc1^2+m2*(l1^2+lc2^2+2*l1*lc2*c2);
d12=I2+m2*(lc2^2+l1*lc2*c2);
d21=d12;
d22=I2+m2*lc2^2;
h=m2*l1*lc2*sin(theta(2));

C=[h*omega(2)*(2*omega(1)+omega(2)); %This represents a sign change.
    h*omega(1)^2];

D=[d11, d12;
    d21, d22];