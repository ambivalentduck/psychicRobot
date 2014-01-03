clc
clear all

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

S=zeros(size(comb.AB,1),1);
ST=S;

N=length(comb.A);

for k=1:length(S)
    S(k)=1/(2*N)*sum((comb.A-comb.AB(k,:)).^2);
    ST(k)=1/N*sum(comb.B.*(comb.AB(k,:)-comb.A));
end

vT=var(comb.AB(:));

names=paramsPopulator('names');
dat=paramsPopulator('burdet');
figure(6)
clf
H=barh(1:20,[S ST 0*S]/vT)
set(H(1),'facecolor',.2*[1 1 1])
set(H(2),'facecolor',.8*[1 1 1])
set(gca,'ytick',1:20)
set(gca,'yticklabel',names(dat(:,4)>0))
xlabel('Fraction of Total Variance')
a1=gca;
a2=axes('Position',get(gca,'position'),'xaxislocation','top','yaxislocation','right','color','none');
set(a1,'ylim',[.25 20.25])
set(a2,'ylim',get(a1,'ylim'))
set(a1,'xlim',[0 1])
%set(a2,'xlim',[0 1])
set(a2,'ytick',[])
set(a1,'box','off')

hold on
for k=1:20
    plot(1000*abs(comb.A-comb.AB(k,:)),k-.5+.15*(rand(1,2000)-.5),'k.','markersize',.0001)
end
set(get(a2,'xlabel'),'string','Change in MUE, mm','verticalalignment','middle')
