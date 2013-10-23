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
    if ty==3
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
    try
        lumps=findLumps(trials(k).t',y,0);
        trials(k).lumps=lumps;
        trials(k).lumpfail=1;
        col=lumps2rgbk(lumps,y);
        for kk=1:1:length(trials(k).t)
            plot(y(kk,1)-y(1,1),y(kk,2)-y(1,2),'.','Color',col(kk,:),'MarkerSize',1)
        end
    catch
        trials(k).lumpfail=0;
        plot(y(1:spacer:end,1)-y(1,1),y(1:spacer:end,2)-y(1,2),'c')
    end
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

