clc
clear all

%% Get data formatted and ready.

if ~exist('OATdat_all.mat','file')
    load all_white_workspace.mat

    white.oat=OAT;
    white.xvaf=xvaf;
    white.yex=yex;

    load baselines.mat
    load OAT_KICK.mat
    load KICK.mat
    load KICKsm.mat

    kick.oat=OAT;

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

    kick.xvaf=xvaf;
    kick.yex=yex;

    save('OATdat_all.mat','kick','white')

    clc
    clear all
end

load('OATdat_all.mat')

%% Calculate MUEs for plotting.

bins=0:.005:.15;
[trash,kick.ref]=getMUE(bins,0*bins,kick.yex);
[trash,white.ref]=getMUE(bins,0*bins,white.yex);

%l1,shoulder xy, sensor bias xy, impedance stuff

dat=paramsPopulator('burdet');
f=find(dat(:,3)>0)

names=paramsPopulator('names');
for k=1:29
    disp([num2str(k),' ',names{f(k)}])
end

whiteVals=zeros(29,9);
kickVals=whiteVals;
for k=1:29
    for kk=1:9
        whiteVals(k,kk)=getMUE(bins,white.ref,white.oat(k,kk).y);
        kickVals(k,kk)=getMUE(bins,kick.ref,kick.oat(k,kk).y);
    end
end

%% Make the plots.

figure(1)
clf
set(0,'defaulttextinterpreter','none')
hold on

oI=[1:2 10:11 16:17 12:13 27:29 18:26];
LRskip=10;
UDskip=-.01;
kickGray=.2;
whiteGray=.8;
prepend='';

xlist=[1:9 9 1];

names=paramsPopulator('latex');

for k=1:20
    mWhite=max(whiteVals(oI(k),:));
    mKick=max(kickVals(oI(k),:));
    if mWhite>mKick
        hf=fill((mod(k-1,4)+1)*LRskip+xlist,(floor((k-1)/4)+1)*UDskip+[whiteVals(oI(k),:) 0 0],kickGray*[1 1 1]);
        hf=fill((mod(k-1,4)+1)*LRskip+xlist,(floor((k-1)/4)+1)*UDskip+[kickVals(oI(k),:) 0 0],whiteGray*[1 1 1]);
    else
        hf=fill((mod(k-1,4)+1)*LRskip+xlist,(floor((k-1)/4)+1)*UDskip+[kickVals(oI(k),:) 0 0],whiteGray*[1 1 1]);
        hf=fill((mod(k-1,4)+1)*LRskip+xlist,(floor((k-1)/4)+1)*UDskip+[whiteVals(oI(k),:) 0 0],kickGray*[1 1 1]);
    end
    %plot((mod(k-1,4)+1)*LRskip+[1 9],(floor((k-1)/4)+1)*UDskip+[0 0],'k','linewidth',2)
    vals=paramsPopulator(f(oI(k)));
    vals=vals(:,f(oI(k)))';
    for kk=[1 5 9]
        th=text((mod(k-1,4)+1)*LRskip+kk,(floor((k-1)/4)+1)*UDskip-.0002,[prepend,num2str(vals(kk),2)],'Horizontalalignment','center','Verticalalignment','top');
    end
    tp=get(th,'Extent');
    s=regexprep(names{f(oI(k))},'Burdet ','');
    s=regexprep(s,', unitless','');
    s=[upper(s(1)) s(2:end)];
    s=[prepend,s];
    text((mod(k-1,4)+1)*LRskip+5,(floor((k-1)/4)+1)*UDskip-.001,s,'fontname','Helvetica','Horizontalalignment','center','Verticalalignment','top')
end

for k=1:5
plot([10 10]-.008,k*UDskip+[.005,0],'k','linewidth',3)
text(10-.008,k*UDskip+.0025,'5 cm MUE','Horizontalalignment','center','Verticalalignment','bottom','rotation',90)
end

ylim([-.055 -.005])
xlim([9 51])
axis on
set(gca,'position',[0 0 1 1],'units','normalized')
axis off
set(gcf,'position',[680, 70. 695, 903])

laprint(gcf,'fig5raw','width',15,'scalefonts','off','factor',1)