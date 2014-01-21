function setGlobals(params)

global l1 l2 m1 m2 lc1 lc2 I1 I2 x0 kp0 kp1 kp kd kpgain kpkdratio reflexratio kpkdreflexratio

v=params;
l1=v(1);
l2=v(2);
lc1=l1*v(3);
lc2=l2*v(4);
mass=v(5);
m1=v(6)*mass;
m2=v(7)*mass;
I1=m1*(v(8)*l1)^2;
I2=m2*(v(9)*l2)^2;

x0=[v(10) v(11)];

kpgain=v(16);
kp0=v(17)*[v(21) v(22); v(23) v(24)];
kp1=v(18)*[v(25) v(26); v(27) v(28)];
kp=v(19)*[v(29) v(30); v(31) v(32)];
kd=v(20)*[v(33) v(34); v(35) v(36)];

kpkdratio=v(37);
reflexratio=v(38);
kpkdreflexratio=v(39);