function out=fJTSKQF(theta)
theta1=theta(1);
theta2=theta(2);
fJ11=- (889.*sin(theta1 + theta2))./2500 - (127.*sin(theta1))./400;
fJ21=(889.*cos(theta1 + theta2))./2500 + (127.*cos(theta1))./400;
fJ12=-(889.*sin(theta1 + theta2))./2500;
fJ22=(889.*cos(theta1 + theta2))./2500;
out=[fJ11,fJ12;fJ21,fJ22];
