% This is a function that takes in an output file and plots error(on the y axis) against
% trial(on the x axis) for baseline, extraction, and the secondary
% baseline(washout) all next to each other for all 8 patients
clc
clear all

% builds one display for all of the folloing data 
figure(1)
% builds the verticle separating lines
for g=1:4
    line([g*96,0],[g*96,0.25])
end
clf
hold on 

%for all 7 patients (1 hasn't been recorded yet)
for num = 1:8
%     changes the number indicie to a string
    truepatientnum=num2str(num+22);
%     concatination with integer string with output and .mat
    output = ['./Data/output',truepatientnum,'.mat']
%     loads the patient data file
    load(output);
%     retrieves the origin and the target
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

% % % %code for the analyzing baselines compared to extraction
%     baseline
%     important(1:96,1)= errors(1:96,1);
%     extract
%     important(97:2*96,1) = errors(1+2*96:3*96,1);
%     washout
%     important(1+2*96:3*96,1) = errors(1+4*96:5*96,1);
    
    xaxis = 1:1:5*96;
    phase3 = xaxis(1+2*96:3*96);
    phase2 = xaxis(97:2*96);
    phase4 = xaxis(3*96+1:4*96);
    phase234 = xaxis(97:4*96);
    subi0=[trials.targcat]~=0;
%     subi1=[trials.targcat]~=1;
    curs0 = subi0(1,1+2*96:3*96);
%     curs1 = subi1(1,1+2*96:3*96);
%     cursx = xaxis(subi(1+2*96:3*96));
%     figure(num)
%     clf
col=rand(3,1);
    plot(xaxis(subi0),errors(subi0,1)+num*.2,'.','color',col)
    dcol=max([col-.3,zeros(3,1)]');
    plot(phase3(curs0),errors(phase3(curs0),2)+num*.2,'x','color',dcol)
    plot(phase3(curs0),mean(errors(phase3(curs0),2))+zeros(48,1)+num*.2,'--','color',dcol)
    plot(phase3(curs0),mean(errors(phase3(curs0),1))+zeros(48,1)+num*.2,'-.','color',col)
    
    errave = (mean(errors(phase2(curs0),1))+zeros(48,1)+num*.2 + mean(errors(phase4(curs0),1))+zeros(48,1)+num*.2)/2;
    plot(phase234(curs0),errave,'--','color','red')
    %axis ([0 5*96 0 .2])
end