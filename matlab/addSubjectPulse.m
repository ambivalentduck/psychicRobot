function addSubjectPulse(k)

subnums=[324 789 300 301 5:8];

name=num2str(subnums(k));
input=load(['../Data/Data_pulse/input',name,'.dat']);
output=load(['../Data/Data_pulse/output',name,'.dat']);
params=getSubjectParams(name);
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

%Because trial 1 starts from an unknown location/where the robot is
%turned on
trials(1).dist=.15;
trials(1).dir=-1;
trials(1).disturbcat=0;

cats=[1 0 2];
dists=[.15 .3];

for c=2:length(trials)
    [trash,ind]=min(abs(norm(trials(c).x(end,:)-trials(c).x(1,:))-dists));
    trials(c).dist=dists(ind);
    trials(c).dir=sign(trials(c).x(end,1)-trials(c).x(1,1));
    trials(c).disturbance=input(c,[4 5 6]);
    v=find(input(c,[4 5 6]));
    if isempty(v)
        trials(c).disturbcat=0;
    else
        trials(c).disturbcat=2*v-(input(c,3+v)>=0);
    end
end
unique([trials.dist])
unique([trials.dir])
unique([trials.disturbcat])
save(['../Data/Data_pulse/pulse',num2str(k),'.mat'],'trials','params')


