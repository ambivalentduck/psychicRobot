clc
clear all
close all

outs=zeros(10,4);

for k=1:10
    subname=num2str(k,'%2.2d');
    fname=['free_exp_',subname,'stroke.mat'];
    
    if ~exist(fname,'file')
        load(['./stroke/QQQ7_s',subname,'_d1.mat']);
        x=DDD(:,1:2);
        filtn=64;
        filtType='loess';
        t=0:.005:.005*(size(DDD,1)-1);
        t=t';
        x=[smooth(x(:,1),filtn,filtType) smooth(x(:,2),filtn,filtType)];
        
        v=[gradient(x(:,1)) gradient(x(:,2))]/.005;
        a=[gradient(v(:,1)) gradient(v(:,2))]/.005;
        save(fname,'t','x','v','a');
    else
        load(fname)
    end
    [outs(k,1), outs(k,2), outs(k,3), outs(k,4)]=boltzmannsubjectplots(t,x,v,a,k)
    suplabel(['Stroke Subject ',subname],'t')
    %set(findall(gcf,'type','text'),'fontSize',24)
    set(gcf,'position',[76 11 1195 925])
    print('-dtiff','-r300',['summary',subname,'stroke.tiff'])
    print('-dpng','-r300',['summary',subname,'stroke.png'])
end

save('stroke.mat','outs')