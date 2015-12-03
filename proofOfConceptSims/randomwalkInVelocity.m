clc
clear all

tstep=.005;
t=0:tstep:100;
v=.1;
x=0;

state=zeros(length(t),3);

W=2; %Might be an outer loop later
hellzno=.05;

%Euler because why not and stochastic
for k=1:length(t)
    if v==0
        v=v+.001*randn; %Floating point zero isn't real and causes nasty headaches
    end
    %r=hellzno+(1-2*hellzno)*rand;
    r=.95*rand;
    invcumulativedensity=-W*log(1-r);
    a=sign(rand-.5)*invcumulativedensity/v;
    
    state(k,:)=[x v a];
    v=v+a*tstep-tstep*(1.6*W)*v;
    x=x+v*tstep;
    
    if x>1
        v=-abs(v);
    elseif x<-1
        v=abs(v);
    end
end

figure(88)
clf
for k=1:3
    subplot(3,3,k)
    plot(t,state(:,k))
    
    subplot(3,3,3+k)
    hist(state(:,k))
end

subplot(3,3,[7])
hist(state(:,2).*state(:,3))

subplot(3,3,8)
hist(log(state(:,2).^2))
