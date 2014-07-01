function [out,rms]=maxperpendicular(x,x0,x1)

% clc
% clear all
% load intern1.mat
% 
% N=15
% x=trials(N).x
% x0=trials(N).orig
% x1=trials(N).targ

x=[x(:,1)-x0(1), x(:,2)-x0(2)];

M=x1-x0;

dists=zeros(size(x,1),1);
for k=1:length(dists)
    dists(k)=norm(x(k,:)-(dot(x(k,:),M)/dot(M,M))*M);
end

[out,ind]=max(dists);

rms=sqrt(mean(dists.^2));