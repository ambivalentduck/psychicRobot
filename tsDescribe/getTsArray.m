function [ts,offsetVec]=getTsArray(x,t,x0,xf,N)

buffer=.1;
offsetVec=linspace(buffer,1-buffer,N+2);
offsetVec=offsetVec(2:end-1);

p=rotateProgressError(x,x0,xf);
onset=find(p(:,1)>=buffer,1,'first');
offsets=zeros(N,1);
for k=1:N
    offsets(k)=find(p(:,1)>=offsetVec(k),1,'first');
end

ts=t(offsets)-t(onset);

function out=rotateProgressError(x,x0,x1)

x=[x(:,1)-x0(1), x(:,2)-x0(2)];

xdiff=x1-x0;

theta=0-atan2(xdiff(2),xdiff(1));
ct=cos(theta);
st=sin(theta);
rotmat=[ct st;-st ct];

out=x*rotmat;
out(:,1)=out(:,1)/norm(xdiff);
