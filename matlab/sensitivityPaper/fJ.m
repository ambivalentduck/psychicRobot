function out=fJ(theta,l1,l2)
theta1=theta(1);
theta2=theta(2);
fJ11=-l1.*sin(theta1)-l2.*sin(theta1+theta2);
fJ21=l1.*cos(theta1)+l2.*cos(theta1+theta2);
fJ12=-l2.*sin(theta1+theta2);
fJ22=l2.*cos(theta1+theta2);
out=[fJ11,fJ12;fJ21,fJ22];
