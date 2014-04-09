clc
clear all

figure(14)
clf
hold on

figure(13)
clf
hold on

for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'U.mat'])

    %In theory, you need to sort by reachcat, which is in reachinfo

    %For now ignore that loop, just sort by subject by baseline, extracted
    %baseline, disturbed hand, extracted from disturbance full and partial,
    %so 5 data per subject.

    clear cleanStruct
    f=find([trialInfo.clean]);
    for t=2:length(f)
        lt=length(trials(f(t)).t);
        lty=length(undisturbed(f(t)).ty);
        inds=(1+lt-lty):lt;
        cleanStruct(t).R2BX=getR2(trials(f(t)).x(inds,:));
        cleanStruct(t).R2BY=getR2(undisturbed(f(t)).y);
    end
    subject(k).R2BX=vertcat(cleanStruct.R2BX);
    subject(k).R2BY=vertcat(cleanStruct.R2BY);

    clear triStruct
    for t=1:length(trials)
        if isempty(trialInfo(t).forceinds)
            continue
        end

        inds=trialInfo(t).forceinds(1)+(0:59);
        onset=find(vecmag(trials(t).v)>.05,1,'first');
        start=max(onset-35,1);
        indsy=max(1,trialInfo(t).forceinds(1)-start)+1+(0:59);
        plot(trials(t).x(inds,1),trials(t).x(inds,2),'k')
        plot(trials(t).y(indsy,1),trials(t).y(indsy,2),'r')
        plot(trials(t).yp(1:length(inds),1),trials(t).yp(1:length(inds),2),'b')

        triStruct(t).R2X=getR2(trials(t).x(inds,:));
        triStruct(t).R2Y=getR2(trials(t).y(indsy,:)); %Alignment issues here
        triStruct(t).R2YP=getR2(trials(t).yp(1:length(inds),:));
    end
    subject(k).R2X=vertcat(triStruct.R2X);
    subject(k).R2Y=vertcat(triStruct.R2Y);
    subject(k).R2YP=vertcat(triStruct.R2YP);

    if k==20000
        blah=subject(k).R2YP;
        [trash,inds]=sort(blah(:,1));
        subject(k).R2YP=subject(k).R2YP(inds(4:end),:);


        %         figure(12)
        %         clf
        %         hold on
        %         plot(trials(t).x(:,1),trials(t).x(:,2),'b',trials(t).y(:,1),trials(t).y(:,2),'r',trials(t).yp(:,1),trials(t).yp(:,2),'g')
        %         axis equal
        %        figure(14)
    end


    CBX=1*ones(size(subject(k).R2BX));
    CBY=2*ones(size(subject(k).R2BY));
    CX=3*ones(size(subject(k).R2X));
    CY=4*ones(size(subject(k).R2Y));
    CYP=5*ones(size(subject(k).R2YP));

    Xlist=[CBX;CBY;CX;CY;CYP];
    Ylist=[subject(k).R2BX;subject(k).R2BY;subject(k).R2X;subject(k).R2Y;subject(k).R2YP];

    xoff=.3/2*((k-1)-1.5);
    figure(14)
    plot(Xlist+xoff,Ylist(:,1),'.','markersize',.0001)
    boxplot(Ylist(:,1),Xlist,'orientation','vertical','notch','on','positions',(1:5)+xoff,'widths',.1,'symbol','r.','labels',{'','','','',''})
end
ylabel('Max Perpendicular Deviation, cm')
set(gca,'xtick',[1:5])
set(gca,'xticklabels',{'Undisturbed','Extracted','Disturbed','Extraction','Extraction$_O$'})
xlim([.7 5.3])
ylim([0 17])
ys=ylim;
plot(2.5*[1 1],ys,'m-')
ylim(ys)
cleanup

set(0,'defaulttextinterpreter','none')
set(gcf,'color',[1 1 1])
set(gcf,'position',[639   176   661   485])
laprint(gcf,'fig4raw','scalefonts','off','asonscreen','on')
