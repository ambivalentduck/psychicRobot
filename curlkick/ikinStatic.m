function theta=ikinStatic(x,l1,l2,x0)

x=x-x0;

theta2=-2*atan(sqrt(((l1+l2)^2-(x(1)^2+x(2)^2))/(x(1)^2+x(2)^2-(l1-l2)^2)));  %+/- at the beginning of this is elbow up/down
theta1=atan2(x(2),x(1))-atan2(l2*sin(theta2),l1+l2*cos(theta2));

theta=[theta1; theta2];
