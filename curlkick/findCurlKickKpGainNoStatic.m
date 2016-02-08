function findCurlKickKpGainNoStatic(k)

load(['../Data/curlkick/curlkick',num2str(k),'.mat'])

%GWR variable controls the width of the gaussian kernel that averages/smooths
GWR=2000;

%% Use undisturbed examples for a reference trajectory for each start/end pair.
dcats=[trials.disturbcat];
reachcat=[trials.targetcat];
rearrangecat=[1 2 0 3];
reachcat=rearrangecat(reachcat)';

nclean=3;
clean=0*dcats;
for kk=1:nclean
    clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
end
clean=(clean==0);

figure(1)
clf
hold on

figure(2)
clf
hold on

toi=0:.01:1;
DOI=2;
for U=DOI
    figure(1)
    f=find(clean'&(reachcat==U));
    rp=randperm(length(f));
    for krp=1:length(rp)
        fk=rp(krp);
        if 1 %krp<=20
            start=curlKickOnsetDetector(trials(f(fk)));
            time=trials(f(fk)).t-trials(f(fk)).t(start);
            x=trials(f(fk)).x;
            pl=sum(sqrt(sum((x(2:end,:)-x(1:end-1,:)).^2,2)));
            if pl>.25 %Misstart and/or wrong target detection
                continue
            end
            catme(U,krp).t=time;
            catme(U,krp).x=[x trials(f(fk)).v trials(f(fk)).a];
        end
        plot(trials(f(fk)).x(:,1),trials(f(fk)).x(:,2),'linewidth',.01)
        plot(trials(f(fk)).x(end,1),trials(f(fk)).x(end,2),'rx')
    end
    catx=vertcat(catme(U,:).x);
    catt=vertcat(catme(U,:).t);
    X=gaussianWeightedRegression(catt,catx,toi,GWR);
    direction(U).X=X;
    plot(X(:,1),X(:,2),'g')
    figure(2)
    plot(catt,vecmag(catx(:,3:4))+.2*(U-1),'b.','markersize',.01)
    plot(toi,vecmag(X(:,3:4))+.2*(U-1),'g')
end
figure(1)
axis equal


%% Get a composite trajectory for each disturbed direction/disturbance pair, 6 in all

figure(3)
clf
hold on

for U=DOI
    for D=2
        toffset=0; %4*(D-1);
        xoffset=.4*(D-1);
        
        figure(3)
        f=find((dcats'==D)&(reachcat==U));
        rp=randperm(length(f));
        for krp=1:length(rp)
            fk=rp(krp);
            if 1 %krp<=20
                start=curlKickOnsetDetector(trials(f(fk)));
                time=trials(f(fk)).t-trials(f(fk)).t(start);
                x=trials(f(fk)).x;
                
                catus(U,D,krp).t=time;
                catus(U,D,krp).x=[x trials(f(fk)).v trials(f(fk)).a trials(f(fk)).f];
            end
            plot(trials(f(fk)).x(:,1)+xoffset,trials(f(fk)).x(:,2),'linewidth',.01)
            plot(trials(f(fk)).x(end,1)+xoffset,trials(f(fk)).x(end,2),'rx')
        end
        catx=vertcat(catus(U,D,:).x);
        catt=vertcat(catus(U,D,:).t);
        X=gaussianWeightedRegression(catt,catx,toi,GWR);
        disturbed(U,D).X=X;
        disturbed(U,D).Y=direction(U).X;
        plot(X(:,1)+xoffset,X(:,2),'g')
        figure(2)
        plot(catt+toffset,vecmag(catx(:,3:4))+.2*(U-1),'r.','markersize',.01)
        plot(toi+toffset,vecmag(X(:,3:4))+.2*(U-1),'c')
    end
end
figure(3)
axis equal
error

X=vertcat(disturbed.X);
Y=vertcat(disturbed.Y);

[l1,l2,x0]=getSubjectParams(k);
set2dGlobals(l1,l2,x0,1); %Setting mass to 1 kg.

[TBMP,TSP,TNP]=cart2model(X,Y);

figure(55)
clf
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on
for kk=5:25
    i=5:(5+kk);
    B=[[TSP(i,1);TSP(i,2)] [TBMP(i,1);TBMP(i,2)]]\[TNP(i,1);TNP(i,2)];
    subplot(2,1,1)
    plot(kk,B(1),'.')
    subplot(2,1,2)
    plot(kk,B(2),'.')
end
%% Relate torques to find shift as well as mass and stiff gains
addpath ../../DeOpt/
params

max([X(:,2);Y(:,2)])
min([X(:,2);Y(:,2)])

[mass,c]=deoptMassStiffShift(TBMP,TSP,TNP,5:20)
error
params.c=c;
params.L1=L1;
params.L2=L2;
params.x0=x0;
params.mass=mass;

save(['../Data/curlkick/curlkick',num2str(k),'W.mat'],'trials','params')
