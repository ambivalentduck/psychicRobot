function out=robotfJ(theta)
theta1=theta(1);
theta2=theta(2);
fJ11=-46251./100000.*sin(theta1)-33521./100000.*sin(theta1+theta2);
fJ21=46251./100000.*cos(theta1)+33521./100000.*cos(theta1+theta2);
fJ12=-33521./100000.*sin(theta1+theta2);
fJ22=33521./100000.*cos(theta1+theta2);
out=[fJ11,fJ12;fJ21,fJ22];
