clc
clear all

name='300'
load(['../Data/',name,'.mat'])

global kpgain

set2dGlobals(params.l1, params.l2, params.origin, params.shoulder, params.mass)

%% Set up vectors for find() later.

early=[trials.early];
late=[trials.late];
white=[trials.white];
disturbed=[trials.updown];

ad=abs(disturbed);
nclean=4;
clean=zeros(size(ad));
for k=1:nclean
    clean=clean+[zeros(1,k-1) ad(1:end-(k-1))];
end
clean=clean==0;

direction=[trials.direction];
dist=[trials.dist];

u=unique([trials(2:end).direction])
u2=unique([trials(2:end).dist])

%% Plot position domain data

spacer=20;
gray=.7;
kpgain=.75;

type=early*pi+late*exp(1)+white; %Use irrational numbers to ease sorting
types=unique(type);
types=types(types~=0);

for F=1:length(types)
    figure(F)
    clf
    for k=1:2
        for kk=1:2
            subplot(2,2,3-k+2*(kk-1))
            hold on
            f=find((direction==u(k))&(dist==u2(kk))&clean);
            for c=1:min(100,length(f))
                plot(trials(f(c)).x(1:spacer:end,1)-trials(f(c)).x(1,1),trials(f(c)).x(1:spacer:end,2)-trials(f(c)).x(1,2),'Color',[gray gray gray])
            end
            plot(0,0,'rx')
            plot([0 cos(u(k))*u2(kk)],[0 0],'m')
            axis equal
        end
    end
end

for k=1:length(type)
    ty=find(types==type(k));
    if isempty(ty)
        continue
    end
    figure(ty)

    f=find(direction(k)==u);
    f2=find(dist(k)==u2);

    subplot(2,2,3-f+2*(f2-1))
    y=extract(trials(k).t',[trials(k).x trials(k).v trials(k).a trials(k).f],@armdynamics_inverted);
    trials(k).extracted=y;
    plot(trials(k).x(1:spacer:end,1)-trials(k).x(1,1),trials(k).x(1:spacer:end,2)-trials(k).x(1,2),'k')
    quiver(trials(k).x(1:spacer:end,1)-trials(k).x(1,1),trials(k).x(1:spacer:end,2)-trials(k).x(1,2),trials(k).f(1:spacer:end,1),trials(k).f(1:spacer:end,2),'b')
    plot(y(:,1)-y(1,1),y(:,2)-y(1,2),'r')
end

figure(1)
suplabel('Kick Early, Left of Rhumb','t')
figure(2)
suplabel('Kick Late, Left of Rhumb','t')
figure(4)
suplabel('Kick Late, Right of Rhumb','t')
figure(5)
suplabel('Kick Early, Right of Rhumb','t')
figure(3)
suplabel('BLWN','t')

%% Plot perpendicular distance against time

%First detect onset of movement, demonstrate validity of detection
figure(1000)
clf
subplot(1,2,1)
hold on
subplot(1,2,2)
hold on

f=find((direction==u(1))&(dist==u2(1))&clean);
for k=1:length(f)
    subplot(1,2,1)
    plot(trials(f(k)).x(:,1)-trials(f(k)).x(1,1),trials(f(k)).v(:,1))

    subplot(1,2,2)
    i=find(vecmag(trials(f(k)).v)>.11,1,'first');
    plot(trials(f(k)).t-trials(f(k)).t(i),trials(f(k)).x(:,1)-trials(f(k)).x(1,1))
end
suplabel('Empirical Determination of Onset of Movement criteria','t')
subplot(1,2,1)
ylabel('Vel_x, m/s')
xlabel('Pos_x, m')
subplot(1,2,2)
ylabel('Pos_x, m')
xlabel('Aligned Time, s')


% Now, slide time around using the above and bin it
binSpan=1; %1 second
binCount=21;
H=(binCount-1)/binSpan;
bins=linspace(0,binSpan,binCount);

%Choose a transform H(t) such that round(H(t)) does the binning.
%H(t)=t*binCount/binSpan
% round((binCount/binSpan)*t)/(binCount/binSpan) -> Binned times

for k=1:length(type)
    i=find(vecmag2(trials(k).v)>=(.11^2),1,'first');
    if isempty(i)
        i=1;
    end
    t=trials(k).t-trials(k).t(i);
    bt=round(H*t)/H;
    for kk=1:length(bins)
        catme(k,kk).vals=trials(k).x(bt==bins(kk),2)';
        if type(k)~=0
            catme(k,kk).exvals=trials(k).extracted(bt==bins(kk),2)';
        end
    end
end

for k=1:2
    for kk=1:2
        f=find((direction==u(k))&(dist==u2(kk))&clean);
        for kkk=1:length(bins)
            catted(1,3-k+2*(kk-1),kkk).vals=[catme(f,kkk).vals];
            catted(1,3-k+2*(kk-1),kkk).mean=mean(catted(1,3-k+2*(kk-1),kkk).vals);
            catted(1,3-k+2*(kk-1),kkk).std=std(catted(1,3-k+2*(kk-1),kkk).vals);
        end

        for T=1:length(types)
            f=find((direction==u(k))&(dist==u2(kk))&(type==types(T)));
            for kkk=1:length(bins)
                catted(T+1,3-k+2*(kk-1),kkk).vals=[catme(f,kkk).vals];
                catted(T+1,3-k+2*(kk-1),kkk).mean=mean(catted(T+1,3-k+2*(kk-1),kkk).vals);
                catted(T+1,3-k+2*(kk-1),kkk).std=std(catted(T+1,3-k+2*(kk-1),kkk).vals);
                catted(T+1,3-k+2*(kk-1),kkk).exvals=[catme(f,kkk).exvals];
                catted(T+1,3-k+2*(kk-1),kkk).exmean=mean(catted(T+1,3-k+2*(kk-1),kkk).exvals);
                catted(T+1,3-k+2*(kk-1),kkk).exstd=std(catted(T+1,3-k+2*(kk-1),kkk).exvals);
            end
        end
    end
end

for F=1:length(types)
    figure(F+10)
    clf
    for k=1:2
        for kk=1:2
            subplot(2,2,3-k+2*(kk-1))
            hold on
            errorbar(bins,[catted(1,3-k+2*(kk-1),:).mean],[catted(1,3-k+2*(kk-1),:).std],'Color',[gray gray gray])
            errorbar(bins+.015,[catted(F+1,3-k+2*(kk-1),:).mean],[catted(F+1,3-k+2*(kk-1),:).std],'k')
            errorbar(bins-.015,[catted(F+1,3-k+2*(kk-1),:).exmean],[catted(F+1,3-k+2*(kk-1),:).exstd],'r')
            xlabel('Movement Time, t')
            ylabel('Pos_x, m')
        end
    end
end

figure(11)
suplabel('Kick Early, Left of Rhumb','t')
figure(12)
suplabel('Kick Late, Left of Rhumb','t')
figure(14)
suplabel('Kick Late, Right of Rhumb','t')
figure(15)
suplabel('Kick Early, Right of Rhumb','t')
figure(13)
suplabel('BLWN','t')
