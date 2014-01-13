clc
clear all

if ~exist('justMUE4fig2.mat','file')
    if ~exist('varSens_all.mat','file')
        load all_white_workspace.mat

        white.yex=yex;
        white.xvaf=xvaf;
        white.A=simA;
        white.B=simB;
        white.AB=simAB;

        t=0:.005:2;
        coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
        tcalc=t;
        tcalc(t>=.7)=.7;
        [x,v,a]=minjerk(coeff,tcalc);
        x=x';
        v=v';
        a=a';

        f=zeros(length(t),2);

        f((t>=.1)&(t<=.15),2)=15;

        xvaf=[x v a f];

        load baselines.mat
        load KICK.mat

        kick.xvaf=xvaf;
        kick.yex=yex;
        kick.A=simA;
        kick.B=simB;
        kick.AB=simAB;

        save('varSens_all.mat','white','kick');
        clc
        clear all
    end

    load varSens_all.mat

    if ~isfield(white,'mueAB')

        bins=0:.005:.15;

        [trash,kick.ref]=getMUE(bins,0*bins,kick.yex);
        [trash,white.ref]=getMUE(bins,0*bins,white.yex);

        white.mueA=zeros(size(white.A));
        white.mueB=zeros(size(white.B));
        white.mueAB=zeros(size(white.AB));
        kick.mueA=zeros(size(kick.A));
        kick.mueB=zeros(size(kick.B));
        kick.mueAB=zeros(size(kick.AB));

        for k=1:length(kick.A)
            white.mueA(k)=getMUE(bins,white.ref,white.A(k).y);
            white.mueB(k)=getMUE(bins,white.ref,white.B(k).y);
            kick.mueA(k)=getMUE(bins,kick.ref,kick.A(k).y);
            kick.mueB(k)=getMUE(bins,kick.ref,kick.B(k).y);
            for kk=1:size(kick.AB,1)
                white.mueAB(kk,k)=getMUE(bins,white.ref,white.AB(kk,k).y);
                kick.mueAB(kk,k)=getMUE(bins,kick.ref,kick.AB(kk,k).y);
            end
        end
        save('varSens_all.mat','white','kick');
    end

    comb.A=[white.mueA kick.mueA];
    comb.B=[white.mueB kick.mueB];
    comb.AB=[white.mueAB kick.mueAB];

    save('justMUE4fig2.mat','comb')
else
    load justMUE4fig2.mat
end

S=zeros(size(comb.AB,1),1);
ST=S;

N=length(comb.A);

for k=1:length(S)
    S(k)=1/(2*N)*sum((comb.A-comb.AB(k,:)).^2);
    ST(k)=1/N*sum(comb.B.*(comb.AB(k,:)-comb.A));
end

vT=var(comb.A(:));

names=paramsPopulator('latex');
dat=paramsPopulator('burdet');
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