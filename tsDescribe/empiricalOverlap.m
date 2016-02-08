clc
clear all

C1=0;
S1=1;

t=-.5:.01:1.5;

td=t-C1; %C1=0
ta=td/S1+.5; %S1=1
K1=(30*ta.^2-60*ta.^3+30*ta.^4)/S1;
K1=K1.*((ta>0)&(ta<1));

S2=1.1;
%C2=.5;

figure(1)
clf
hold on

k=0;
%for C2=0:.01:(C1+S1/2+S2/2)
for C2=.01:.01:5
    td=t-C2;
    ta=td/S2+.5;
    K2=(30*ta.^2-60*ta.^3+30*ta.^4)/S2;
    K2=K2.*((ta>0)&(ta<1));
    
    %extra=2*K1.*K2.*(t>(C2-S2/2))&(t<(C1+S1/2));
    extra=sqrt(2*K1.*K2);
    E=trapz(t,extra);
    
    k=k+1;
    %plot3(k*ones(size(t)),t,extra)
    %plot(C2,E,'k.')
    plot(C2,E,'k.')
end

%plot(t,extra)