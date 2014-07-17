function addSubject(k)

% if k is in a character form then that's what the name will be
if isa(k,'char')
    name=k;
    c=1;
else
    name=['output',num2str(k)];
    c=2;
end

% % % % % THE INPUT FILE IS WORKING WITH THE OLD INPUT!!!!
input=load('./Data/input.dat');
output=load(['./Data/',name,'.dat']);

% if - elseif statement that tells what getSubjectParams to look for
% if c==1 (the first case) the user has put in a string. Therefore the
% parameters for that particular subject must be under that specific
% name
% if c==2 (the second case) the user has put in a number. Therefore the
% parameters must be under a name that is an integer (which is then converted
% into a string)
if c==1
    params=getSubjectParams(name);
elseif c==2
    params=getSubjectParams(num2str(k));
end

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
y=y(inds,:);

if size(output,2)>=14
    cursor=output(:,13:14);
    cursor=cursor(inds,:);
else
    cursor=y;
end


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

% % % Why is the following code that is commented here?

% if isa(k,'char')
%     name=k;
% else
%     name=['intern',num2str(k)];
% end

save(['./Data/',name,'.mat'],'trials','params')
disp(name)
