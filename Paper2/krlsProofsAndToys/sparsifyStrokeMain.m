clc
clear all

global xdot t

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

%Set global
xdot=v;

%find a point with high speed and low curvature (pray that means it's
%pretty pure)
[val,ind]=max(K.*vm);

%Initialize a "center" lump at that point and zero on either side a decent
%guess at lump width away
ini=[v(ind,:)'*.7/1.875;t(ind);.7;0;0;0;t(ind)-.35;.7;0;0;0;t(ind)+.35;.7];

%Perform the center optimization and plot it
P=fminunc(@sparseupMJPfit,ini,optimoptions(@fminunc,'GradObj','on','TolX',1E-9));

%In each direction, do 2 unknown, 1 known until the magnitude comes back
%tiny.


r=reshape(P,5,length(P)/5);
w=r(1:3,:);
tc=r(4,:);
ts=abs(r(5,:));

uppers=tc+ts/2;
lowers=tc-ts/2;

inds=find((t>lowers(1))&(t<uppers(1)));
tinds=t(inds);

kern=zeros(length(tinds),3);
summed=zeros(length(tinds),3);
for k=1:3
    kern(:,k)=slmj5op(tinds,tc(k),ts(k));
    summed=summed+(w(:,k)*kern(:,k)')';
end

figure(2)
clf
plot(tinds,vecmag(summed),'r',tinds,vecmag(xdot(inds,:)),'k',tinds,vecmag(summed-xdot(inds,:)),'g')
