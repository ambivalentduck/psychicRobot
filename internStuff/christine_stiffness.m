clc
clear all

phases=zeros(96*5,1);
for k=1:5
    phases(96*(k-1)+1:96*k)=k;
end

load './Data/christine.mat'

global K

set2dGlobals(params.l1, params.l2, params.shoulder,params.mass) %Make sure mass actually in kg

F=find(phases==3);

for k=1:length(F)
    xvaf=[trials(F(k)).x trials(F(k)).v trials(F(k)).a trials(F(k)).f];
    qt=xvaf2arm(xvaf);

    x0=xvaf(1,1:2); %orig;
    x1=trials(F(k)).targ;
    x=[xvaf(:,1)-x0(1), xvaf(:,2)-x0(2)];
    M=x1-x0;

    ft=find((trials(F(k)).t-trials(F(k)).t(1))<.5,1,'last');
    
    xRhumb=xvaf(1:ft,1:2);
    qRhumb=xRhumb;
    xdat(k).t=qt(1:ft,7:8);
    
    for kk=1:ft
        xRhumb(kk,:)=xRhumb(kk,:)-(x(kk,1:2)-(dot(x(kk,:),M)/dot(M,M))*M);
        qRhumb(kk,:)=ikin(xRhumb(kk,:));
        
        %[D,C]=computeDC(qt(kk,1:2),qt(kk,3:4));
        %xdat(k).t(kk,:)=qt(kk,7:8); %+(D*qt(kk,7:8)'-C)'-;
    end
    
    xdat(k).q=qRhumb(:,1:2)-qt(1:ft,1:2);

end

q=vertcat(xdat.q);
t=vertcat(xdat.t);

K=q\t

figure(10)
clf
hold on

for k=1:length(F)
    xvaf=[trials(F(k)).x trials(F(k)).v trials(F(k)).a trials(F(k)).f];
    trials(F(k)).yfixed=extract(trials(F(k)).t,xvaf,@armdynamics_inverted);
    plot(xvaf(:,1),xvaf(:,2),'k',trials(F(k)).y(:,1),trials(F(k)).y(:,2),'r',trials(F(k)).yfixed(:,1),trials(F(k)).yfixed(:,2),'b')
end
axis equal



figure(1)
clf
plot(q(:,1),t(:,1),'.','markersize',1e-6)
%axis equal

figure(2)
clf
plot(q(:,2),t(:,2),'.','markersize',1e-6)
%axis equal

cleanup