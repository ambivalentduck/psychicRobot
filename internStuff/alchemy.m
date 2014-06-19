clc
clear all
load intern3.mat

for N=1:length(trials)
    x=trials(N).x;
    y=trials(N).y;
    x0=trials(N).orig;
    x1=trials(N).targ;
    
    [xout, xind]=maxperpendicular(x,x0,x1);
    [yout, yind]=maxperpendicular(y,x0,x1);
    
    dists(N,1)=xout;
    dists(N,2)=yout;
end

disp(dists)

close all

figure(1)
clf
plot(dists(:,1),dists(:,2),'.')

f=find((dists(:,1)<.5)&(dists(:,2)<.5));

[dists(f,1) 0*dists(f,1)+1]\dists(f,2)

phases=zeros(96*5,1);
for k=1:5
    phases(96*(k-1)+1:96*k)=k;
end

[p,blah,stats]=anova1(dists(:,1)-dists(:,2),phases);
multcompare(stats)

X=dists(:,1);
X(96*2+1:96*3)=dists(96*2+1:96*3,2);

augphases=[phases; 6*ones(96,1)];

[p,blah,stats]=anova1([X;dists(96*2+1:96*3,1)],augphases);
multcompare(stats)