clc
clear all

%% Step 0: Prove this is a RKHS

if 0 %If you've run this once, there's no compelling reason to do it again.
    fh=figure(500)
    clf
    t=0:.005:1;
    W=.4;
    exemp=.5;
    v=slmj5op(t,exemp,W);
    for tk=1:length(t)
        clf
        hold on
        plot(t,v,'k')
        plot(t,slmj5op(t(tk),t,W),'r')
        plot(exemp,slmj5op(exemp,t(tk),W),'b.','markersize',5)
        F(tk)=getframe;
    end
    movie(fh,F,1,50)
end


%% Step 1: Make a sum of mj polynomial submovements.

tfinal=.7;
firstspan=tfinal/2;

fspansc=5/6;
tc=[firstspan/2 tfinal/2 tfinal-firstspan/2]
ts=[firstspan fspansc*tfinal firstspan];

cols=.8*rand(3,length(tc));

ti=tc-ts/2;
tf=tc+ts/2;

x0=[1 1 zeros(1,4)];
sMag=.015;
dx=[sMag sMag;.15 0;sMag -sMag];

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


%% Step 2: See how kernel RLS turns out.

%Initialize the only kernel parameter: its span.
%Fundamental tradeoff is that you want these about as wide as you can
%manage without ever breaking the fundamental rule: no more than
%2-coverage.
ks=.25; %So at any given time, we're carrying around ks/.005=30 kernels. Wide span risks overfitting, but probably gets the direction right

%H in time is weird because you already perfectly know FUTURE input samples, but not output samples.
%Oddly this implies that H never changes because of symmetry.
bound=floor(ks/(2*t_step))*t_step;
H=slmj5op(0,-bound:t_step:bound,ks);
H_old=H(1:end-1);
H_new=H(end);
H0=slmj5op(0,0,ks);
lH=length(H);

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
W=zeros(lH-1+length(t),2);
%a is a subset of W.

d=yg(:,3:4)+.05*randn(length(t),2);
e_mag=zeros(length(t),1);

for iter=1:length(t)
    e=d(iter,:)-H'*W(iter:iter+lH-1,:);
    e_mag(iter)=norm(e);
    %2/sum(H) is...bizarre but very clearly optimal.
    W(iter:iter+lH-1,:)=W(iter:iter+lH-1,:)+sH2*H_reg.*H*e;
end

%% Step 3: In theory we just learned something.
t_expanded=t(1)-t_step*(lH-1)/2:t_step:t(end)+t_step*(lH-1)/2;

recon=zeros(length(t),2);
figure(1)
[summed,kerns2]=supMJP([1 1 0 0 0 0],W,t_expanded-ks/2,t_expanded+ks/2,t);
subplot(2,1,1)
hold on
plot(summed(:,1),summed(:,2),'r.')
title('Position')


axis equal
subplot(2,1,2)
hold on
plot(t,sqrt(sum(d.^2,2)),'k')
plot(t,sqrt(sum(summed(:,3:4).^2,2)),'r.')
xlabel('Time, s')
ylabel('Speed, m/s')
title('Speed')
for k=1:length(tc)
    plot(t,norm(dx(k,:))*kerns(:,k),'color',cols(:,k))
end
legend('Ground Truth','Noisy Speed Supplied to Algorithm','Reconstruction','Submovement 1','Submovement 2','Submovement 3')


% figure(99)
% clf
% plot(t,sqrt(sum(d.^2,2))-sqrt(sum(summed(:,3:4).^2,2)))


subset=(lH-1)/2+(1:length(t));

% figure(2)
% clf
% subplot(2,1,1)
% plot(t,e_mag)
% ylabel('Innovation (Prediction Error), m/s')
% xlabel('Time, s')
% subplot(2,1,2)
% plot(t,atan2(W(subset,2),W(subset,1)))
%
% figure(3)
% plot(t,gradient(atan2(W(subset,2),W(subset,1))))

%% Step 4: Intent with zero noise, verify tautology

global M B
M=.41;
B=2.3;

F=([0; 6]*((t>.1)&(t<.150)))';
K=([15; 15]*ones(size(t)))';

xva=forward(t,yg,F,K+2);
figure(5)
clf
subplot(2,1,1)
hold on
plot(yg(:,1),yg(:,2),'g-','linewidth',5)
plot(xva(:,1),xva(:,2),'k')
arrow(xva(:,1:2),xva(:,1:2)+F/100,[.5 .5 .5],.3);
axis equal

taut=reverse(t,xva,F,K);
plot(taut(:,1),taut(:,2),'m')

%% Can we back out stiffness?

for k=1:3
    out{k}=reverse(t,xva,F,K+k);
end

ks=.25;
lambda=1.1;

bound=floor(ks/(2*t_step))*t_step;
H=slmj5op(0,-bound:t_step:bound,ks);
H_old=H(1:end-1);
H_new=H(end);
H0=slmj5op(0,0,ks);
lH=length(H);
sH2=3.5/sum(H);

H_reg=H;
H_reg(1:(lH-1)/2)=0; %The trailing edge should not cause learning
H_reg=H_reg/sum(H_reg);
H_reg((lH-1)/2+1:end)=H_reg(end:-1:(lH-1)/2+1);

%W is larger than t because we're going to have kernels that touch the start and end of movement.
W=zeros(lH-1+length(t),2);

e_mag=zeros(length(t),1);

dX=(out{3}(:,3:4)-out{2}(:,3:4));
d=out{1}(:,3:4)+16*dX;
Ke=10+zeros(length(t),2);

P=1E3*eye(2);

for iter=2:length(t)
    %Calculate error
    e=d(iter,:)-(H'*W(iter:iter+lH-1,:)+Ke(iter-1,:).*dX(iter,1:2));
    e_mag(iter)=norm(e);

    %Update the kernels
    W(iter:iter+lH-1,:)=W(iter:iter+lH-1,:)+sH2*H_reg.*H*e;

    %Update the parameters
    x=dX(iter,1:2)';
    k=(lambda*P*x)/(1+lambda*x'*P*x);
    Ke(iter,:)=Ke(iter-1,:)+k'.*e; %Diag-only lets us .*, usually a matrix.
    P=lambda*P-lambda*k*x'*P;

    %LMS-style update fails, but leave the code here to show examples if need be
    %Ke(iter,:)=Ke(iter-1,:)+e_mag(iter)*100*lambda*dX(iter,:).*e;
end

[summed,kerns2]=supMJP([1 1 0 0 0 0],W,t_expanded-ks/2,t_expanded+ks/2,t);
plot(summed(:,1),summed(:,2),'r.')

subplot(2,1,2)
hold on
plot(t,yg(:,4),'g','linewidth',5)
plot(t,xva(:,4),'k')
%plot(t,out{1}(:,4),'m')
plot(t,10*dX(:,2),'color',[.5 .5 .5])
%plot(t,d(:,2)+14*dX(:,2),'m')
plot(t,summed(:,4),'r.')
%plot(t,summed(:,4)+Ke(:,2).*dX(:,2),'r--')
xlabel('Time,s')
ylabel('Velocity in Cartesian Y direction, m/s')