clc
clear all

if ~exist('justMUE4fig2.mat','file')
    loadme=1:8;

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
            if ~mod(c,100)
                plot(simA(c).y(:,1),simA(c).y(:,2),'g')
            end
            mueA(c)=getMUE(bins,ref,simA(c).y);
            mueB(c)=getMUE(bins,ref,simB(c).y);
            for cc=1:20
                mueAB(c,cc)=getMUE(bins,ref,simAB(cc,c).y);
            end
        end
        justMUE(k).mueA=1000*mueA;
        justMUE(k).mueB=1000*mueB;
        justMUE(k).mueAB=1000*mueAB;
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
        justMUE(k).rawS(kk,:)=(abs(justMUE(k).mueA-justMUE(k).mueAB(:,kk))/2)';
        justMUE(k).rawST(kk,:)=sqrt(abs(justMUE(k).mueB.*(justMUE(k).mueAB(:,kk)-justMUE(k).mueA)))';
    end
end

allS=[justMUE.rawS];
allST=[justMUE.rawST];

lrs=length(allS);
rp=randperm(lrs);
plotS=allS(:,rp(1:1000));
plotST=allST(:,rp(1:1000));

allA=vertcat(justMUE.mueA);
allB=vertcat(justMUE.mueB);
allAB=vertcat(justMUE.mueAB);

N=8000;
for k=1:20
    ST(k)=1/(2*N)*sum((allA-allAB(:,k)).^2);
    S(k)=1/N*sum(allB.*(allAB(:,k)-allA));
end
v=[S' ST']
vall=var(allAB(:))
sub_white=allAB(1:4000,:);
v_white=var(sub_white(:))
sub_pulse=allAB(4001:8000,:);
v_pulse=var(sub_pulse(:))

v=v/vall
totalSensitivity=sum(sum(v))
%v=v/totalSensitivity;
[vals,order]=sort(v(:,2));
%order=1:length(order)

height=8;
width=6;
lmargin=8;
bmargin=1;
figmargin=.4;

figure(1)
clf
hold on
set(gcf,'color',[1 1 1])
%set(gcf,'units','centimeters')
%set(gcf,'position',[4,8,figmargin+lmargin+width,figmargin+bmargin+height])
%set(gcf,'units','normalized')

spacer=.1;
NPLOT=1000;
scatter=.15*(rand(1,NPLOT)-.5);
rawS=[justMUE.rawS];
rawST=[justMUE.rawST];
R=randperm(size(rawS,2));
oNPLOT=ones(1,NPLOT);
MSIZE=2;
for k=1:20
    plot(rawS(order(k),R(1:NPLOT)),(k+spacer)*oNPLOT+scatter,'.','markersize',MSIZE,'color',[.5 0 .5])
    plot(rawST(order(k),R(1:NPLOT)),(k-spacer)*oNPLOT+scatter,'.','markersize',MSIZE,'color',[0 .8 .8])
end

names=paramsPopulator('latex');
dat=paramsPopulator('burdet');
names=names(dat(:,4)>0);

set(gca,'xcolor',[1 1 1])
set(gca,'ytick',1:20)
set(gca,'yticklabel',names(order))
set(gca,'xtick',[])

TIGHT=.5;
ylim([TIGHT 21-TIGHT])

%title('Change in MUE')
text(21,-2,'Change in MUE')

labN=.05;
plot([0 1]+labN,TIGHT*[1 1],'k','linewidth',4)
text(labN+.5,TIGHT,'1 mm','Horizontalalignment','center','verticalalignment','top')

margin=.05;

set(gcf,'position',[250 300 625 400])
set(gca,'position',[.5 margin 1-margin 1-margin],'units','normalized')

set(0,'defaulttextinterpreter','none')
set(gcf,'color',[1 1 1])

laprint(gcf,'../figures/fig2raw','scalefonts','off','asonscreen','on')

%print ../figures/fig2 -depsc

for k=1:20
    disp([num2str(v(k,1),'%0.2f'),' ',num2str(v(k,2),'%0.2f'),' ',names{k}])
end

