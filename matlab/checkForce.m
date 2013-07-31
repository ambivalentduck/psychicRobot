clc
clear all
close all

load ../Data/output11.dat
load ../Data/input11.dat
load ../Data/11.mat

set2dGlobals(params.l1,params.l2,params.origin,params.shoulder,params.mass)

a=unique(input11(:,2));
val=sqrt(sum(output11(:,[5 6]).^2,2))+sqrt(sum(output11(:,[7 8]).^2,2));
inTarget=zeros(size(val));
for k=[1 2 4]
    inTarget=inTarget|abs(output11(:,3)-a(k))<.0125;
end
f=find((val<.1)&inTarget);

figure(1)
subplot(2,2,1)
hold on
plot(output11(f,7),output11(f,9),'.')
plot(output11(f,8),output11(f,10),'r.')
ylabel('Force Sensor, N')
xlabel('Handle Acceleration, m/s^2')
legend('X','Y')
subplot(2,2,2)
filt=ones(1,1)/1;
hold on
plot(output11(f,1),output11(f,9),'b.',output11(f,1),output11(f,10),'r.',output11(f,1),filter(filt,1,output11(f,9)),'c-',output11(f,1),filter(filt,1,output11(f,10)),'m-','LineWidth',3)

kick=find(sum(abs(input11(max(1,output11(f,1)-1),[4 5])),2));
plot(output11(f(kick),1),output11(f(kick),10),'kx')
xlabel('Trial')
ylabel('Force Sensor Offset, N')
legend('X','Y','X Trend','Y Trend','Kick Came Just Before')

o=params.origin;
op=ikin(o);
e=[filter(filt,1,output11(f,9))';filter(filt,1,output11(f,10))'];
kp=SnMKpGainStruct(1.5);
kp=kp{1};
diff=f;
normf=f;
o=o';
for k=1:length(f)
    diff(k)=norm(o-fkin(kp\e(:,k)+op));
    normf(k)=norm(e(:,k));
end
[counts,bins]=hist(diff);
subplot(2,2,3)
plot(output11(f,1),diff*100,'.')
xlabel('Trial')
ylabel('Estimated extraction offset due to force sensor drift, cm')

subplot(2,2,4)
plot(normf,diff*100,'.')
xlabel('Force Sensor Drift *Magnitude*, N')
ylabel('Extraction Error *Magnitude*, cm')

figure(2)
N=5;
o=filter(ones(N,1)/N,1,output11(:,10));
%o=filter([1 1]/2,[1 1]/2,output11(:,10));
o=o(f);
xo=xcorr(o);
xo=xo/sum(xo);
x=xcorr(output11(f,10));
x=x./sum(x);
filt=ones(50,1)/50;
r=rand(size(output11(f,10)));
x2=xcorr(filter(filt,[1],r));
x2=x2./sum(x2);
x3=xcorr(filter(filt,[1],r)+100*rand(size(output11(f,10))));
x3=x3./sum(x3);
i=-floor(length(x)/2):floor(length(x)/2);
plot(i,x,i,xo,'x',i,x2,i,x3)
title('Autocorrelation')
legend('Real data','Filtered Real','Simulated, no noise','Simulated+Noise')

figure(3)
plot(filter(filt,1,[zeros(50,1);1; zeros(50,1)]))

figure(4)
a=unique(output11(f,1));
off=0;
offx=0;
trialoff=NaN*ones(length(a),1);
trialoffx=NaN*ones(length(a),1);
for k=1:length(a)
    i=find(output11(f,1)==a(k));
    trialoff(a(k))=mean(output11(f(i),10));
    trialoffx(a(k))=mean(output11(f(i),9));
end
for k=1:length(trialoff)
    if isnan(trialoff(k))
        trialoff(k)=off;
        trialoffx(k)=offx;
    else
        off=trialoff(k);
        offx=trialoffx(k);
    end
end
for k=2:length(trialoff)
    trialoffd(k)=trialoff(k)-trialoff(k-1);
    trialoffdx(k)=trialoffx(k)-trialoffx(k-1);
end


off=0;
trialoff2=NaN*ones(length(a),1);
for k=1:length(a)
    i=find(output11(f,1)==a(k));
    trialoff2(a(k))=mean(o(i));
end
for k=1:length(trialoff2)
    if isnan(trialoff2(k))
        trialoff2(k)=off;
    else
        off=trialoff2(k);
    end
end
for k=2:length(trialoff2)
    trialoffd2(k)=trialoff2(k)-trialoff2(k-1);
end

x1=xcorr(trialoffd);
x1=x1/sum(abs(x1));
x2=xcorr(trialoffd2);
x2=x2/sum(abs(x2));
x=xcorr(output11(f,10));
x=x./sum(x);
i=-floor(length(x1)/2):floor(length(x1)/2);
plot(i,x1,i,x2,'.')
title('Autocorrelation')
legend('No prefilter','Short prefilter')

figure(5)
clf
hold on
plot(output11(f,3),output11(f,4),'.')
a=unique(input11(:,2));
t=0:.01:2*pi;
circlex=cos(t);
circley=sin(t);
SCALE=.025/2;
for k=[1 2 4]
    plot(SCALE*circlex+a(k),SCALE*circley+.5,'r')
end
axis equal

figure(6)
o=params.origin;
op=ikin(o);
e=[trialoffdx;trialoffd];
kp=SnMKpGainStruct(1.5);
kp=kp{1};
diff=f;
normf=f;
o=o';
diff=trialoff;
for k=1:length(trialoff)
    diff(k)=norm(o-fkin(kp\e(:,k)+op));
end
plot(1:length(trialoff),diff*100,'.')
counts2=hist(diff(diff~=0));
xlabel('Trial')
ylabel('Steady-State Error')
title('Corrected Error by Trial')

figure(7)
bar(bins*100,[counts/sum(counts);counts2/sum(counts2)]'); %apples and oranges
ylabel('Frequency')
xlabel('Steady-state Error, cm')
legend('Without correction','Subtracting off the most recent clump')
