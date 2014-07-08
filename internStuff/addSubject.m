function addSubject(k)

if isa(k,'char')
    name=k;
else
    name=['output',num2str(k)];
end

input=load('./Data/old_input.dat');
output=load(['./Data/',name,'.dat']);
params=getSubjectParams(name);

%outStream << trial TAB now-zero TAB position.X() TAB position.Y() TAB velocity.X() TAB velocity.Y() TAB accel.X() TAB accel.Y() TAB force.X() TAB force.Y() TAB desposition.X() TAB desposition.Y() TAB xpcTime << endl;

traw=output(:,2);
xvafraw=output(:,3:10);
trial=output(:,1);
lastTrial=max(trial)-1;
f=find(trial==lastTrial,1,'last');
samprate=1/mean(gradient(traw))
if samprate<300 %our filtering method apparently has nonlinear algorithmic complexity and goes from seconds to hours+.
    inds=1:f;
else
    inds=1:5:f;
end
trials=generalAddSubject(name,traw(inds),xvafraw(inds,:),trial(inds),params);

y=output(:,11:12);
cursor=output(:,13:14);
y=y(inds,:);
trial=trial(inds);

%Because trial 1 starts from an unknown location/where the robot is
%turned on

cats=[1 0 2];
dists=[.15 .3];

[a,b,c]=unique(sum(input(:,2:3),2));
remap=[1 0 2 3];
targs=c;

targlocs=input(b,2:3);


for c=1:length(trials)
    trials(c).disturbance=input(c,6);
    if c==1
        trials(c).orig=trials(c).x(1,:);
    else
        trials(c).orig=targlocs(targs(c-1),:);
    end
    trials(c).targ=targlocs(targs(c),:);
    trials(c).targcat=remap(targs(c));
    trials(c).y=y(trial==c,:);
    trials(c).cursor=cursor(trial==c,:);
end

if isa(k,'char')
    name=k;
else
    name=['intern',num2str(k)];
end

save(['./Data/',name,'.mat'],'trials','params')
