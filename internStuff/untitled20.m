clc
clear all

nS=1;
nE=2;
lT=5*96;

error.x=zeros(lT,nE,nS);
error.y=zeros(lT,nE,nS);
error.f=zeros(lT,nE,nS);
error.t=zeros(lT,1,nS);

phases=zeros(96*5,1);
for k=1:5
    phases(96*(k-1)+1:96*k)=k;
end

%% Load and process

S=1

if ~exist('Cerror.mat','file')|1
%     load(['./Data/christine21.mat']);
      load(['./Data/output22.mat']);

    for N=1:length(trials)
        x=trials(N).x;
        y=trials(N).y;
        f=trials(N).f;
        x0=trials(N).orig;
        x1=trials(N).targ;

        onset=find(vecmag(trials(N).v)>.05,1,'first');
        start=max(onset-35,1);
        ft=find((trials(N).t-trials(N).t(start))<inf,1,'last');
        
        start=start+15;
        
        [xmp, xrms]=maxperpendicular(x(start:ft,:),x0,x1);
        [ymp, yrms]=maxperpendicular(y(start:ft,:),x0,x1);
        [fmp, frms]=maxperpendicular(f(start:ft,:),x0,x1);

        error.x(N,:,S)=[xmp,xrms];
        error.y(N,:,S)=[ymp,yrms];
        error.f(N,:,S)=[fmp,frms];
        error.t(N,1,S)=trials(N).t(end)-trials(N).t(1);
    end

    save('Cerror.mat','error')
else
    load('Cerror.mat');
end

close all

nV=1:lT;

mX=mean(error.x,3);
mY=mean(error.y,3);
mT=mean(error.t,3);

mXF=mean(error.f./error.x,3);
mYF=mean(error.f./error.y,3);

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

[p,blah,stats]=anova1([mX(:,1); mY(phases==3,1)],[phases; 0*phases(phases==3)+3.5]);
multcompare(stats)

[p,blah,stats]=anova1([mXF(:,2); mYF(phases==3,2)],[phases; 0*phases(phases==3)+3.5]);
multcompare(stats)


return