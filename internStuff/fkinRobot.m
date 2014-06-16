function [x,x1]=fkinRobot(theta)

l1=.46251;
l2=.33521;

x=[l1*cos(theta(1))+l2*cos(theta(1)+theta(2));
   l1*sin(theta(1))+l2*sin(theta(1)+theta(2))];

if nargout>1
    x1=[l1*cos(theta(1));
        l1*sin(theta(1))];
end

