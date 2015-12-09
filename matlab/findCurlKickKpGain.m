function findCurlKickKpGain(k)

load(['../Data/curlkick/curlkick',num2str(k),'.mat'])

global kpgain massgain

kpgain=1;
massgain=1;

params

GWR=2000;

%% Categorize by start/end pair

%Irrelevent as the target unique identifies the origin of outward reaches
%now

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
for U=1:3
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
    X=gaussianWeightedRegression(catt,catx(:,1:4),toi,GWR);
    plot(X(:,1),X(:,2),'g')
    figure(2)
    plot(catt,vecmag(catx(:,3:4))+.2*(U-1),'b.','markersize',.01)
    plot(toi,vecmag(X(:,3:4))+.2*(U-1),'g')
end
figure(1)
axis equal


%% Use a helper function to get [x v a], [xd vd ad], and F into torques

%This code is super straight-forward: error and error velocity vs typical
%error and typical error velocity. F=k1(Da+C)+k2(Ze)
% inertia and coriolis are linear in body mass
% k2 is essentially a stiffness adjustment
figure(3)
clf
hold on
error
fitlower=.0/.005;
fitupper=.6/.005;
indspan=fitlower:fitupper;
for U=1:length(urc)
    f=find((dcats==5)&(reachcat'==(urc(U))));
    catx=vertcat(catme(U,:).x);
    catt=vertcat(catme(U,:).t);
    for fk=1:length(f)
        kk=f(fk);
        
        start=onsetDetector(trials(kk));
        time=trials(kk).t(start+indspan);
        time=time-time(1);
        
        %Find typical reach-perpendicular p,v,a.
        
        Y=gaussianWeightedRegression(catt,catx(:,1:6),time,GWR);
        Y(:,2)=Y(:,2)-(Y(1,2)-trials(f(fk)).x(start+indspan(1),2)); %If the "shape" holds, this corrects for starting bias.
        
        X=[trials(kk).x(start+indspan,:) trials(kk).v(start+indspan,:) trials(kk).a(start+indspan,:) trials(kk).f(start+indspan,:)];
        
        plot(time,X(:,1)+yoffset*U,'b')
        plot(time,Y(:,1)+yoffset*U,'r')
        
        [storeme(U,fk).TMP,storeme(U,fk).TSP,storeme(U,fk).TNP]=cart2model(X,Y);
        storeme(U,fk).time=time;
        storeme(U,fk).X=X;
        storeme(U,fk).Y=Y;
        storeme(U,fk).starts=starts(kk);
        storeme(U,fk).ends=ends(kk);
    end
end
axis equal

%% Relate torques to find shift as well as mass and stiff gains

scat=[storeme.starts];
subfitlower=.25/.005;
subfitupper=.45/.005;
for kk=1:3
    f=find(scat);
    TMP=[storeme(f).TMP];
    TSP=[storeme(f).TSP];
    TNP=[storeme(f).TNP];
    [massgains(kk),kpgains(kk),lag]=deoptMassStiffShift(TMP,TSP,TNP,subfitlower:subfitupper)
end


%% Shift force trace relative to position trace due to sensor lead.
if 0
    lag=round(lag);
    f_all=vertcat(trials.f);
    
    for T=3:length(trials)-2 %Quietly fail to fix the first two and last two reaches because they don't matter and are undisturbed
        trials(T).f=f_all(trials(T).rawinds+lag,:); %#ok<*AGROW>
    end
end

save(['../Data/Data_pulse/pulse',num2str(k),'W.mat'],'massgains','kpgains','trialInfo','means','trials','params','storeme')
