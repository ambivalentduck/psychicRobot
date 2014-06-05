clc
clear all

name='301'
load(['../../Data/',name,'.mat'])


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

u=unique([trials(2:end).direction]);
u2=unique([trials(2:end).dist]);

f=find(clean&(direction==u(1))&(dist==u2(1)));

for k=1:length(f)
    catme(k).x=[trials(f(k)).x(:,1)-trials(f(k)).origin(1) trials(f(k)).x(:,2)-trials(f(k)).origin(2)]';
    catme(k).t=trials(f(k)).target(1);
end
x=[catme.x]';
figure(1)
clf
hold on
plot(x(:,1),x(:,2),'b.')
plot([0 .15],[0 0],'rx')
axis equal

for k=1:length(f)
    mueReal(k)=getMUE(bins,ref2,[trials(f(k)).x(:,1)-trials(f(k)).origin(1) trials(f(k)).x(:,2)-trials(f(k)).origin(2)]);
end

var(mueReal)

load simSobol.mat
load baselines.mat

mueSim=zeros(size(simmedAB));

iyex=interp1(yex(:,1)+.0000001*rand(size(yex(:,1))),yex(:,2),bins);

for k=1:size(simmedAB,1)
    for kk=1:size(simmedAB,2)
        i=interp1(simmedAB(k,kk).y(:,1)+.0000001*rand(size((simmedAB(k,kk).y(:,1)))),simmedAB(k,kk).y(:,2),bins);
        i=i-iyex;
        i=i(~isnan(i));
        mueSim(k,kk)=mean(abs(i))*1000;
    end
end

var(mueSim(:))

for k=1:10000
    R=randperm(length(mueSim(:)));
    v73(k)=var(mueSim(R(1:73)));
end
figure(2)
hist(v73)

load sobolSM.mat

mueSim=zeros(size(saltelliAB));

iyex=interp1(yexsm(:,1)+.0000001*rand(size((yexsm(:,1)))),yex(:,2),bins);

for k=1:size(saltelliAB,1)
    for kk=1:size(saltelliAB,2)
        i=interp1(saltelliAB(k,kk).y(:,1)+.0000001*rand(size((saltelliAB(k,kk).y(:,1)))),saltelliAB(k,kk).y(:,2),bins);
        i=i-iyex;
        i=i(~isnan(i));
        mueSim(k,kk)=mean(abs(i))*1000;
    end
end

var(mueSim(:))

