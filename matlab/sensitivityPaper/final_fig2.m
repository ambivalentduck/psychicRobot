clc
clear all

if ~exist('justMUE4fig2.mat','file')
    whites=6:10;
    whites=whites(randperm(5));
    loadme=[whites(1:4) 11:14];

    bins=0:.005:.15;

    setGlobals(paramsPopulator)

    for k=1:length(loadme)
        load(['BATCH',num2str(loadme(k)),'.mat'])
        yex=extract(t,xvaf,'reflex');

        [trash,ref]=getMUE(bins,0*bins,yex);

        mueA=zeros(1000,1);
        mueB=zeros(1000,1);
        mueAB=zeros(1000,20);

        for c=1:1000
            mueA(c)=getMUE(bins,ref,simA(c).y);
            mueB(c)=getMUE(bins,ref,simB(c).y);
            for cc=1:20
                mueAB(c,cc)=getMUE(bins,ref,simAB(cc,c).y);
            end
        end
        justMUE(k).mueA=mueA;
        justMUE(k).mueB=mueB;
        justMUE(k).mueAB=mueAB;
    end

    save('justMUE4fig2.mat','justMUE')
else
    load justMUE4fig2.mat
end

for k=1:8
    justMUE(k).rawS=zeros(20,1000);
    justMUE(k).rawST=zeros(20,1000);
    justMUE(k).total=std(justMUE(k).mueAB(:));
    for kk=1:20
        justMUE(k).rawS(kk,:)=(1000*abs(justMUE(k).mueA-justMUE(k).mueAB(:,kk))/2)';
        justMUE(k).rawST(kk,:)=1000*sqrt(abs(justMUE(k).mueB.*(justMUE(k).mueAB(:,kk)-justMUE(k).mueA)))';
    end
end

allS=[justMUE.rawS];
allST=[justMUE.rawST];

lrs=length(allS);
rp=randperm(lrs);
plotS=allS(:,rp(1:1000));
plotST=allST(:,rp(1:1000));

% for k=1:length(S)
%     S(k)=1/(2*N)*sum((comb.A-comb.AB(k,:)).^2);
%     ST(k)=1/N*sum(comb.B.*(comb.AB(k,:)-comb.A));
% end

figure(6)
clf
hold on
spacer=.1;
NPLOT=1000;
scatter=.15*(rand(1,NPLOT)-.5);
rawS=[justMUE.rawS];
rawST=[justMUE.rawST];
R=randperm(size(rawS,2));
oNPLOT=ones(1,NPLOT);
MSIZE=2;
for k=1:20
    plot(rawS(k,R(1:NPLOT)),(k+spacer)*oNPLOT+scatter,'.','markersize',MSIZE,'color',[.5 0 .5])
    plot(rawST(k,R(1:NPLOT)),(k-spacer)*oNPLOT+scatter,'.','markersize',MSIZE,'color',[0 .5 .5])
end

names=paramsPopulator('latex');
dat=paramsPopulator('burdet');

set(gca,'ytick',1:20)
set(gca,'yticklabel',names(dat(:,4)>0))

ylim([0 21])

xlabel('Change in MUE, mm')

set(0,'defaulttextinterpreter','none')
laprint(gcf,'fig2raw','scalefonts','off','asonscreen','on')

return



figure(6)
clf
set(gcf,'position',[375   261   629   540])
hold on
spacer=.17;
width=.1;
H1=barh((1:20)+1.5*spacer,abs(S)/vT,'barwidth',width)
H2=barh((1:20)-.5*spacer,abs(ST)/vT,'barwidth',width)
set(H1,'facecolor',.2*[1 1 1])
set(H2,'facecolor',.5*[1 1 1])
set(gca,'ytick',1:20)
set(gca,'yticklabel',names(dat(:,4)>0))
a1=gca;
p1=[.5 .1 .45 .9]; %get(a1,'position');
realYlim=[0 20.5];
a2spacer=2;
h=a2spacer/(realYlim(2)-realYlim(1)+a2spacer);
p2=p1+[0 h 0 -h]*p1(4);
a2=axes('Position',p2,'xaxislocation','bottom','yaxislocation','left','color','none');
set(a1,'ylim',realYlim-[a2spacer 0])
set(a2,'ylim',realYlim)
set(a1,'xlim',[0 1])
%set(a2,'xlim',[0 1])
set(a2,'ytick',[])
set(a1,'box','off')
%set(a1,'ticklength',[0 0])

hold on
for k=1:20
    plot(1000*abs(comb.A-comb.AB(k,:))/2,k+.5*spacer+width*(rand(1,2000)-.5),'.','color',.2*[1 1 1],'markersize',.0001)
    plot(1000*sqrt(abs(comb.B.*(comb.AB(k,:)-comb.A))),k-1.5*spacer+width*(rand(1,2000)-.5),'.','color',.5*[1 1 1],'markersize',.0001)
end
upXLh=annotation('textbox',[.494 .006 .456 .0666])
set(upXLh,'string','Bar Fraction of Total Variance (${\sim}1 {\mu}m$)','edgecolor','none','horizontalalignment','left')
%set(upXLh,'string','Bar Fraction of Total Variance (${\sim}1 {\mu}m$)','horizontalalignment','center')
lXLh=annotation('textbox',[.5435 .072 .456 .0666])
set(lXLh,'string','Change in MUE, mm','edgecolor','none','horizontalalignment','left')
% set(a2xl,'string','Change in MUE, mm','verticalalignment','middle')
% set(a2xl,'position',get(a2xl,'position')+[0 3 0])6
set(gcf,'color','w')
xt=get(a2,'xtick')
set(a2,'xtick',xt(2:end))
set(a1,'position',p1)

set(0,'defaulttextinterpreter','none')
laprint(gcf,'fig2raw','scalefonts','off','asonscreen','on')