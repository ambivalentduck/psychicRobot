function [D,C]=computeDCStatic(theta,omega,l1,l2,mass)

m1=.028*mass;
m2=.022*mass;
lc1=.436*l1;
lc2=.682*l2;
I1=m1*(.322*l1)^2;
I2=m2*(.468*l2)^2;

%Compute alpha to torque relationship, eq. 7.87 in Spong's Robot Control
%and Modeling: pg 262

c2=cos(theta(2));
d11=I1+I2+m1*lc1^2+m2*(l1^2+lc2^2+2*l1*lc2*c2);
d12=I2+m2*(lc2^2+l1*lc2*c2);
d21=d12;
d22=I2+m2*lc2^2;

h=m2*l1*lc2*sin(theta(2)); %sign change. should be -? - reflects coordinate difference. + reflects burdet 2006
C=[h*omega(2)*(2*omega(1)+omega(2));
    h*omega(1)^2];

D=[d11, d12;
    d21, d22];