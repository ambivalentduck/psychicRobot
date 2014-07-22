function dx=toyInvDyn(t,x)

global M B measuredVals measuredTime

interped=twoNearestNeighbor(measuredVals,measuredTime,t);
xh=interped(1:2)';
vh=interped(3:4)';
ah=interped(5:6)';
F=interped(7:8)';
K=interped(9:10)';

xe=x(1:2);
ve=x(3:4);

dx=[ve; (M*ah+B*(vh-ve)+K.*(xh-xe)-F)/M];