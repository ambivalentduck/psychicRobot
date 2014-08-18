function dx=toyDyn(t,x)

global M B measuredVals measuredTime

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
xe=interped(1:2)';
ve=interped(3:4)';
ae=interped(5:6)';
F=interped(7:8)';
K=interped(9:10)';

xh=x(1:2);
vh=x(3:4);

dx=[vh; (M*ae+B*(ve-vh)+K.*(xe-xh)+F)/M];