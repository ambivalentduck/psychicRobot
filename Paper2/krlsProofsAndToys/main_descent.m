clc
clear all
close all


tfinal=.7;
firstspan=tfinal/2;

tc=[firstspan/2 tfinal/2 tfinal-firstspan/2]
ts=[firstspan 5*tfinal/6 firstspan];

cols=.8*rand(3,length(tc));

ti=tc-ts/2;
tf=tc+ts/2;

x0=[1 1 zeros(1,4)];
dx=[.02 .02;.1 0;.02 -.02];

t_step=.005;
t=0:t_step:1;

[yg,kerns]=supMJP(x0,dx,ti,tf,t);

figure(1)
clf
subplot(2,1,1)
plot(yg(:,1),yg(:,2),'g-','linewidth',5)
axis equal

subplot(2,1,2)
hold on
plot(t,sqrt(sum(yg(:,3:4).^2,2)),'g-','linewidth',5)

%% Part 2 is a perfect expansion because why not?
ks=.25; %So at any given time, we're carrying around ks/.005=30 kernels. Wide span risks overfitting, but probably gets the direction right
lambda=0.01; %No reason to think we should punish large velocities, but let's preempt spikes to infinity.

%H in time is weird because you already perfectly know FUTURE input samples, but not output samples.
%Oddly this implies that H never changes because of symmetry.
bound=floor(ks/(2*t_step))*t_step;
H=slmj5op(0,-bound:t_step:bound,ks);
H_old=H(1:end-1);
H_new=H(end);
H0=slmj5op(0,0,ks);
lH=length(H);

Hdiv=H/sum(H);

sH2=3.5/sum(H);
H_reg=H;
H_reg(1:(lH-1)/2)=0; %The trailing edge should not cause learning
H_reg=H_reg/sum(H_reg);
H_reg((lH-1)/2+1:end)=H_reg(end:-1:(lH-1)/2+1);

%In the "degenerate" case, we're already in motion but don't have learned
%coefficients. Ignore this case because we're never forced to deal with
%this in practice. Presume we start from rest.
Q=zeros(lH); %(lambda+H0)*eye(lH);

%W is larger than t because we're going to have kernels that touch the start and end of movement.
W=zeros(length(t),2);

%Because c*0 != 0, weights outside t0+ks/2 < t < tf-ks/2 are uniformly
%zero.
inds=1:length(t);

for k=1:length(t)-50
    i=k:k+50;
    W(k+25,:)=[dot(Hdiv,yg(i,3)) dot(Hdiv,yg(i,4))]/5.7;
end

[rg,kerns]=supMJP(x0,W,t-ks/2,t+ks/2,t);

subplot(2,1,1)
hold on
plot(rg(:,1),rg(:,2),'r-')

subplot(2,1,2)
plot(t,sqrt(sum(rg(:,3:4).^2,2)),'r-')

