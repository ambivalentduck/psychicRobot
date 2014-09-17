clc
clear all

%Is there a way to treat decomposition like a puzzle? Any time v and a are
%both low, you can break the problem into subproblems. Where else can you
%do that?

load outputsS02D1R2_org_session.mat

close all

figure(1)
clf

T=3;

t=Robot_arm(T).time;
x=Robot_arm(T).smooth_pos;
v=Robot_arm(T).vel;
a=[gradient(v(:,1),t) gradient(v(:,2),t) gradient(v(:,3),t)];
vm=sqrt(sum(v.^2,2));
K=t;
for k=1:size(x,1)
    K(k)=norm(cross(v(k,:),a(k,:)))/(vm(k).^3);
end

K=smooth(K);
K=-log(K);
K=K-min(K);
K=K/max(abs(K));


subplot(1,2,1)
plot3(x(:,1),x(:,2),x(:,3))
axis equal

subplot(1,2,2)
plot(t,v(:,1),t,v(:,2),t,v(:,3),t,vm,'k',t,K,'k.',t,K.*vm,'r.')

%Finding points of low curvature and high velocity is equivalent to finding
%the corners of a puzzle?

%If you think you have a center and a direction worth starting from, in
%theory you just need a span.

%In theory, you can make this episodic via optimization 3-deep w,tc,ts and
%cost for making things more than 2-deep. Few free parameters.

