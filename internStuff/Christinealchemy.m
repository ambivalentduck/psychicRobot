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

phases=zeros(96*5,1);
for k=1:5
    phases(96*(k-1)+1:96*k)=k;
end

if ~exist('error.mat','file')
    for S=1:8
        load(['./Data/intern',num2str(S+1),'.mat']);

        for N=1:length(trials)
            x=trials(N).x;
            y=trials(N).y;
%           trials with yfix
            yf = trials(N).yfix;
            f=trials(N).f;
            x0=trials(N).orig;
            x1=trials(N).targ;
            

            dist2lt=.01^2;

            if phases(N) == 3
                dist2=(y(:,1)-x1(:,1)).^2 + (y(:,2)-x1(:,2)).^2;
            else
                dist2=(x(:,1)-x1(:,1)).^2 + (x(:,2)-x1(:,2)).^2;
            end
            lastind=find(dist2 < dist2lt,1,'first');
            

            [xmp, xrms]=maxperpendicular(x(1:lastind,:),x0,x1);
            [ymp, yrms]=maxperpendicular(y(1:lastind,:),x0,x1);
%           maxperpen for yfix
            [ymp, yrms]=maxperpendicular(yf(1:lastind,:),x0,x1);
            [fmp, frms]=maxperpendicular(f(1:lastind,:),x0,x1);

            error.x(N,:,S)=[xmp,xrms];
            error.y(N,:,S)=[ymp,yrms];
%           error for yfix
            error.yf(N,:,S)=[ymp,yrms];
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
% mean for yfix
myf=mean(error.yf,3);
mT=mean(error.t,3);

mXF=mean(error.f./error.x,3);
mYF=mean(error.f./error.y,3);
% error divide for yfix
myfF=mean(error.f./error.yf,3);
mYF=mean(error.f./error.y,3);

figure(1)
clf
hold on
plot(nV,mX(:,1),'k')
plot(nV,mY(:,1),'r')
% nv vs myf plot
plot(nV,myf(:,1),'blue')
plot(nV,mXF(:,2),'ko')
plot(nV,mYF(:,2),'r.')
% nv vs myfF
plot(nV,myfF(:,2),'blue')
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

%modified ttest with myfF for third phase 
[h,p]=ttest(mXF(96+1:2*96,2)-mYF(96+1:2*96,2),mXF(2*96+1:3*96,2)-mYF(2*96+1:3*96,2),myfF(2*96+1:3*96,2)-myfF(2*96+1:3*96,2))


[p,blah,stats]=anova1([mX(:,1); mY(phases==3,1)],[phases; 0*phases(phases==3)+3.5]);
multcompare(stats)

[p,blah,stats]=anova1([mXF(:,2); mYF(phases==3,2)],[phases; 0*phases(phases==3)+3.5]);
multcompare(stats)

%stats for myf?
[p,blah,stats]=anova1([mXF(:,2); mYF(phases==3,2)],[phases; 0*phases(phases==3)+3.6]);
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