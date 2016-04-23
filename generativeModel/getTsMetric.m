function [ts, rw]=getTsMetric(x,v,a,t,x0,xf,varargin)

if nargin<7
    buffer=.1;
else
    buffer=varargin{1};
end

ts=zeros(length(buffer),1);
rw=ts;

p=rotateProgressError(x,x0,xf);

for k=1:length(buffer)
    onset=find(p(:,1)>buffer(k),1,'first');
    offset=find(p(:,1)>(1-buffer(k)),1,'first');
    ts(k)=t(offset)-t(onset);
    rw(k)=trapz(t(onset:offset),abs(dot(v(onset:offset,:)',a(onset:offset,:)')'));
end

function out=rotateProgressError(x,x0,x1)

x=[x(:,1)-x0(1), x(:,2)-x0(2)];

xdiff=x1-x0;

theta=0-atan2(xdiff(2),xdiff(1));
ct=cos(theta);
st=sin(theta);
rotmat=[ct st;-st ct];

out=x*rotmat;
out(:,1)=out(:,1)/norm(xdiff);
