function dx=integrateTorque(t,x)

global MBKTvarying measuredTime

interped=twoNearestNeighbor(MBKTvarying,measuredTime,t);
M=interped(1); 
B=interped(2);
K=interped(3);
T=interped(4);

polyC=conv([M B K],[M B K]);

dx=[x(2:4); (T-(polyC(2)*x(4)+polyC(3)*x(3)+polyC(4)*x(2)+polyC(5)*x(1)))/polyC(1)];