clc
clear all

nS=8;
nE=2;
lT=5*96;

error.x=zeros(lT,nE,nS);
error.y=zeros(lT,nE,nS);
error.f=zeros(lT,nE,nS);
error.t=zeros(lT,1,nS);

%% Load and process

if ~exist('error.mat','file')
    for S=1:8
        load(['./Data/intern',num2str(S+1),'.mat']);

        for N=1:length(trials)
            x=trials(N).x;
            y=trials(N).y;
            f=trials(N).f;
            x0=trials(N).orig;
            x1=trials(N).targ;

            [xmp, xrms]=maxperpendicular(x,x0,x1);
            [ymp, yrms]=maxperpendicular(y,x0,x1);
            [fmp, frms]=maxperpendicular(f,x0,x1);

            error.x(N,:,S)=[xmp,xrms];
            error.y(N,:,S)=[ymp,yrms];
            error.f(N,:,S)=[fmp,frms];
            error.t(N,1,S)=trials(N).t(end)-trials(N).t(1);
        end
    end
    save('error.mat','error')
else
    load('error.mat');
end


%% Plot
close all

nV=1:lT;

mX=mean(error.x,3);
mY=mean(error.y,3);
mT=mean(error.t,3);

mXF=mean(error.x./error.f,3);
mYF=mean(error.y./error.f,3);

figure(1)
clf
hold on
plot(nV,mX(:,1),'k')
plot(nV,mY(:,1),'r')
plot(nV,mXF(:,2),'ko')
plot(nV,mYF(:,2),'r.')
uLim=.1;
plot(2*96+[0 0],[0 uLim],'m')
plot(3*96+[0 0],[0 uLim],'m')
ylim([0 uLim])

figure(2)
clf
hold on
plot(nV,mXF(:,2)-mYF(:,2),'g.')
lLim=-.02;
uLim=.02;
plot(2*96+[0 0],[lLim uLim],'m')
plot(3*96+[0 0],[lLim uLim],'m')
ylim([lLim uLim])

[h,p]=ttest(mXF(96+1:2*96,2)-mYF(96+1:2*96,2),mXF(2*96+1:3*96,2)-mYF(2*96+1:3*96,2))

phases=zeros(96*5,1);
for k=1:5
    phases(96*(k-1)+1:96*k)=k;
end

[p,blah,stats]=anova1(mX(:,1),phases);
multcompare(stats)


return

plot(dists(:,1),dists(:,2),'.')

f=find((dists(:,1)<.5)&(dists(:,2)<.5));

[dists(f,1) 0*dists(f,1)+1]\dists(f,2)



[p,blah,stats]=anova1(dists(:,1)-dists(:,2),phases);
multcompare(stats)

X=dists(:,1);
X(96*2+1:96*3)=dists(96*2+1:96*3,2);

augphases=[phases; 6*ones(96,1)];

[p,blah,stats]=anova1([X;dists(96*2+1:96*3,1)],augphases);
multcompare(stats)