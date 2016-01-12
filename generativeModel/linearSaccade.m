function [x,v]=linearSaccade(varargin)

if nargin==0
    x=1;
    return
end

t0=varargin{1};
tf=varargin{2};
t=varargin{3};

x=zeros(size(t));
v=zeros(size(t));

inds=((t>=t0)&(t<=tf));
x(inds)=(t(inds)-t0)/(tf-t0);
v(inds)=1/(tf-t0);
