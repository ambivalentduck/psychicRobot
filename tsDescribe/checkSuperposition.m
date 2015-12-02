function checkSuperposition

%t must be 1 by n
t=0:.01:3;
tc1=1;
tc2=2;
T1=2;
T2=2;

%L should be D by 1
L1=[1;0]/sqrt(2);
L2=[-1;0]/sqrt(2);

[v1,E1]=compKern(t,tc1,T1,L1);
[v2,E2]=compKern(t,tc2,T2,L2);

figure(1)
clf
hold on
x=cumtrapz(v1+v2,2)'*.01;
plot(t,x(:,1),'b')

z=cumtrapz(sqrt(E1+E2),2)'*.01;
plot(t,z,'r')

end

function [v,E]=compKern(t,tc,T,L)
T2=T/2;
td=t-tc;
ta=td/T+.5;
logind=(td>=-T2)&(td<=T2);
ta=ta(logind);
v=zeros(length(L),length(t));
v(:,logind)=L*(30*ta.^2-60*ta.^3+30*ta.^4)/T;
E=dot(v,v);
end
