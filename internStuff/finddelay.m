clc
clear all

S=26 %why not
if S<20
    prefix='intern';
else
    prefix='output';
end

load(['./Data/',prefix,num2str(S),'.mat'])

a=vertcat(trials(3:479).a);
f=vertcat(trials(3:479).f);

l=size(a,1);

figure(42)
clf
hold on

for blah=3:12
ind=floor(blah*l/16):floor((blah+1)*l/16);
pairs=zeros(31,2);
for k=-30:30
    mdl=fitlm(a(ind,1),f(ind+k,1));
    pairs(31+k,:)=[k*.005 mdl.Rsquared.Ordinary];
end
[v,i]=max(pairs(:,2));
plot(1000*pairs(:,1),pairs(:,2))
plot(1000*pairs(i,1),pairs(i,2),'rx')
end
xlabel('Force Sensor Lag, ms')
ylabel('R^2, regression of acceleration and force')

return 
data=load(['./Data/',prefix,num2str(S),'.dat']);

x=data(:,3:4);
v=[gradient(x(:,1)) gradient(x(:,2))]/.001;
a=[gradient(v(:,1)) gradient(v(:,2))]/.001;
f=data(:,9:10);

l=size(a,1);

for blah=3:12
ind=floor(blah*l/16):floor((blah+1)*l/16);
pairs=zeros(31,2);
for k=-30:30
    mdl=fitlm(a(ind,1),f(ind+k,1));
    pairs(31+k,:)=[k*.005 mdl.Rsquared.Ordinary];
end
[v,i]=max(pairs(:,2));
plot(1000*pairs(:,1),pairs(:,2),'k')
plot(1000*pairs(i,1),pairs(i,2),'gx')
end

