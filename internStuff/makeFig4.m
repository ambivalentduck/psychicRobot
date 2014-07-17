% This is a function that takes in an output file and plots error(on the y axis) against
% trial(on the x axis) for baseline, extraction, and the secondary
% baseline(washout) all next to each other for all 8 patients
clc
clear all
close all

load('./Data/errors.mat')

allT=1:5*96;
phases=zeros(5*96,6);
for k=[1 2 4 5]
    phases(96*(k-1)+1:k*96,k)=1;
end
phases=phases(:,[1 2 3 6 4 5]);
phases(2*96+1:2.5*96,3)=1;
phases(2.5*96+1:3*96,4)=1;

phases=phases==1;

% builds one display for all of the folloing data

units={'Max Error, cm','Stiffness, N/m','Error Span, cm','Time to complete reach, s'};
spacer=[.2 450 .14];

SUBS=8;
colSpan=.5;
subWidth=colSpan/(1.5*SUBS);

hand(:,:,2)=force(:,:,2)./hand(:,:,2);
cursor(:,:,2)=force(:,:,2)./cursor(:,:,2);

% The figure that will display all three plots on one figure


map=[1 4 2 3];

exportme=zeros(5*96,6);

for S=1:SUBS
    inds=5*96*(S-1)+1:5*96*S;
%     exportme(inds,1)=cursor(:,S,1);
%     exportme(inds,2)=cursor(:,S,3);
%     exportme(inds,3)=reachT(:,S);
%     exportme(inds,4)=hand(:,S,2);
    exportme(inds,1)=log10(cursor(:,S,1));
    exportme(inds,2)=log10(cursor(:,S,3));
    exportme(inds,3)=log10(reachT(:,S));
    exportme(inds,4)=log10(hand(:,S,2));
    exportme(inds,5)=phases*(1:6)';
    exportme(inds,6)=S+inds*0;
end
f=find(isnan(exportme(:,4)));
exportme(f,4)=exportme(f-1,4);

fid=fopen('dostats.csv','w');
fprintf(fid,'MaxE,SpanE,ReachT,Stiff,Phase,Sub\n');
fprintf(fid,'%g,%g,%g,%g,%g,%g\n',exportme');
fclose(fid);

for k=1:4
    figure(k)
    [p,table,stats]=anovan(exportme(:,k),exportme(:,[5 6]),'display','off')
    multcompare(stats)
end

figure(10)
clf
set(gcf,'color','w')

for N=1:4
    subplot(1,4,map(N));
    hold on
    
    for P=1:6
    switch N
        case {1 3}
            y=log10(100*cursor(phases(:,P),:,N));
        case 2
            y=log10(hand(phases(:,P),:,N));
        case 4
            y=log10(reachT(phases(:,P),:));
    end
    y=y(~isnan(y));
    m=mean(y);
    st=std(y)/sqrt(length(y));
    
    ALPHA=.3;
    fill([m-st m m+st m], [P P+colSpan/2 P P-colSpan/2],'k','edgealpha',ALPHA,'facealpha',ALPHA)
    
    end
    
    for S = 1:SUBS
        
        col=[.7*rand(1,3)];
        dcol=1-col;
        
        xoff=colSpan*((S-1)/(SUBS-1))-colSpan/2;
        
        for P=1:6
            switch N
                case {1 3}
                    y=log10(100*cursor(phases(:,P),S,N));
                case 2
                    y=log10(hand(phases(:,P),S,N));
                case 4
                    y=log10(reachT(phases(:,P),S));
            end
            sz=sum(phases(:,P));
            plot(y,ones(sz,1)*P+xoff+subWidth*(rand(sz,1)-.5),'.','markersize',.00001,'color',col)
        end
    end
    if N == 2
        %xlabel('Trials')
    end
    switch N
                case {1 3}
                    y=log10(100*cursor(:,:,N));
                case 2
                    y=log10(hand(:,:,N));
                case 4
                    y=log10(reachT(:,S));
    end
    y=y(~isnan(y));
    xlim([quantile(y(:),.01) quantile(y(:),.99)])
    xlabel(units{N})
    if N~=1
        set(gca,'ytick',[])
    end
end
drawnow
%set(gcf,'position',[1821,33,360,859])
subplot(1,4,1)
set(gca,'ytick',1:6)
set(gca,'yticklabel',{'Baseline','Forces','Forces + Intent','Forces + Intent','Forces','Washout'})
set(gcf,'renderer','opengl')