function [maxpd,rms,straightness]=maxperpendicular(x,x0,x1)

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
MdM=dot(M,M);
Mp=[M(2) -M(1)];
Mp=Mp/norm(Mp);

dists=zeros(size(x,1),1);
for k=1:length(dists)
    V=x(k,:)-(dot(x(k,:),M)/MdM)*M;
    dists(k)=norm(V);
    perp(k)=dot(V,Mp);
end

[maxpd,ind]=max(dists);
rms=sqrt(mean(dists.^2));
straightness=max(perp)-min(perp);