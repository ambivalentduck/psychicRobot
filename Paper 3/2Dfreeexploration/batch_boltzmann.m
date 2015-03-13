clc
clear all
close all

sublets={'JL','JT','LP','MF','YM'};

figure(4000)
clf
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on

for k=1:length(sublets)
    fname=['free_exp_',sublets{k},'.mat'];
    
    if ~exist(fname,'file')
        mat=load(['./healthy/h_1fch_ps1_',sublets{k},'.dat']);
        
        %find the first free exploration blocks
        fstop=find((mat(1:end-1,10)==2)&mat(2:end,10)~=2);
        fstart=find((mat(1:end-1,10)~=2)&mat(2:end,10)==2)+1;
        figure(100)
        clf
        hold on
        inds=1:size(mat,1);
        plot(inds,mat(:,10))
        plot(inds(fstop),mat(fstop,10),'rx')
        plot(inds(fstart),mat(fstart,10),'go')
        tinds=fstart(1):fstop(2); %This particular subject
        
        t=0:.005:(mat(tinds(end),11)-mat(tinds(1),11));
        x=mat(tinds,[1 2]);
        
        filtn=64;
        filtType='loess';
        x=[smooth(x(:,1),filtn,filtType) smooth(x(:,2),filtn,filtType)];
        
        v=[gradient(x(:,1)) gradient(x(:,2))]/.005;
        a=[gradient(v(:,1)) gradient(v(:,2))]/.005;
        save(fname,'t','x','v','a');
    else
        load(fname)
    end

    f=find((vecmag(v)>.2)|(vecmag(a)>1));
    t=t(f);
    x=x(f,:);
    v=v(f,:);
    a=a(f,:);
    
    figure(k)
    clf
    [outs(k),P,E]=boltzmannsubjectplots(t,x,v,a)

    suplabel(['Healthy Subject ',sublets{k}],'t')
    set(findall(gcf,'type','text'),'fontname','Times','fontSize',12)
    set(gcf,'position',[76 11 1195 925])
    print('-dtiff','-r300',['summary',sublets{k},'.tiff'])
    print('-dpng','-r300',['summary',sublets{k},'.png'])
    
    figure(4000)
    subplot(2,1,1)
    plot(P(:,2),P(:,1),'b')
    subplot(2,1,2)
    plot(E(:,2),E(:,1),'b')
end

save('healthy.mat','outs')

