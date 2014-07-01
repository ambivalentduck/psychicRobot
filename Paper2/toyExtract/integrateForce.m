function dx=integrateForce(t,x)

global M B measuredTime measuredVals

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
F=interped(7)';
K=interped(9)';

polyC=conv([M B K],[M B K]);

dx=[x(2:4); (F-(polyC(2)*x(4)+polyC(3)*x(3)+polyC(4)*x(2)+polyC(5)*x(1)))/polyC(1)];