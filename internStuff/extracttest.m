% this chunk plots the extraction trials and displays the extraction phase
% (phase 3) with every piece of data that is targcat ~= 0
load('./Data/output23.mat')


% loads up the error data in to a vector called errors 
for N=1:length(trials)

x0=trials(N).orig;
x1=trials(N).targ;
x = trials(N).x;

onset=find(vecmag(trials(N).v)>.05,1,'first');
start=max(onset-35,1);
ft=find((trials(N).t-trials(N).t(start))<inf,1,'last');
[errors(N,1)]=maxperpendicular(x(start:ft,:),x0,x1);
errors(N,2)=maxperpendicular(trials(N).cursor(start:ft,:),x0,x1);
end

% setting paramerters for the plot step
xaxis = 1:1:5*96;
phase3 = xaxis(1+2*96:3*96);
subi0=[trials.targcat]~=0;
%     subi1=[trials.targcat]~=1;
curs0 = subi0(1,1+2*96:3*96);
%     curs1 = subi1(1,1+2*96:3*96);
%     cursx = xaxis(subi(1+2*96:3*96));

% plot step
figure(1)
clf
hold on
plot(x,y,'x','color','r')
plot(phase3,errors(phase3),'o')