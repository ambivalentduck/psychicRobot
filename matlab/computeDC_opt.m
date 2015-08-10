function [D,C]=computeDC_opt(theta, omega,l1,paramvec)

lc1=paramvec(1);
lc2=paramvec(2);
m1=paramvec(3);
m2=paramvec(4);
I1=paramvec(5);
I2=paramvec(6);

%Compute alpha to torque relationship, eq. 7.87 in Spong's Robot Control
%and Modeling: pg 262

c2=cos(theta(2));
d11=I1+I2+m1*lc1^2+m2*(l1^2+lc2^2+2*l1*lc2*c2);
d12=I2+m2*(lc2^2+l1*lc2*c2);
d21=d12;
d22=I2+m2*lc2^2;

h=-m2*l1*lc2*sin(theta(2)); %sign change. should be -? - reflects coordinate difference. + reflects burdet 2006
C=[h*omega(2)*(2*omega(1)+omega(2));
    h*omega(1)^2];

D=[d11, d12;
    d21, d22];