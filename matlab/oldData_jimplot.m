clc
clear all

submax=8;
backgray=.6;

figure(667)
clf

starts={};
ends={};
groupCell={};

nlumps=zeros(45,submax);
SPAN=.8;
cols=colorScheme(10);

for k=1:submax
    load(['../Data/curlkick',num2str(k),'g.mat']);
    targetcat=[trials.targetcat];
    disturbcat=[trials.disturbcat];

    h=hist(targetcat,unique(targetcat));
    [trash,center]=max(h);

    targlocs=zeros(max(targetcat),2);
    for c=1:max(targetcat)
        f=find(targetcat==c);
        x=zeros(length(f),2);
        for cc=1:length(f)
            x(cc,:)=trials(f(cc)).x(end,:);
        end
        targlocs(c,:)=mean(x);
    end
    targdists=vecmag([targlocs(:,1)-targlocs(center,1) targlocs(:,2)-targlocs(center,2)]);

    f=find((targetcat~=center)&(disturbcat~=0));
    nl=[trials(f).nlumps];
    ip=[trials(f).isproblem];

    for c=1:5
        f2=find((nl>=c)&~ip);
        starts{k,c}=zeros(3,length(f2));
        ends{k,c}=zeros(3,length(f2));
        groupCell{k,c}=c*ones(1,length(f2));
        for cc=1:length(f2)
            kk=f(f2(cc));
            inds=trials(kk).lumps(c).inds;
            rat=targdists(targetcat(kk));
            starts{k,c}(:,cc)=[trials(kk).cumdist(inds(1))/rat trials(kk).yurot(inds(1),1)/rat trials(kk).tcat(inds(1))]';
            ends{k,c}(:,cc)=[trials(kk).cumdist(inds(end))/rat trials(kk).yurot(inds(end),1)/rat trials(kk).tcat(inds(end))]';
        end
    end


    f=find((targetcat~=center)&(disturbcat~=0));
    for c=1:length(f)
        kk=f(c);
        if ~trials(kk).isproblem
            for L=1:length(trials(kk).lumps)
                inds=trials(kk).lumps(L).inds;
                rat=targdists(targetcat(kk));
                starts_=[trials(kk).cumdist(inds(1))/rat trials(kk).yurot(inds(1),1)/rat trials(kk).tcat(inds(1))];
                ends_=[trials(kk).cumdist(inds(end))/rat trials(kk).yurot(inds(end),1)/rat trials(kk).tcat(inds(end))];

                for SP=1:3
                    subplot(5,3,3*SP-2:3*SP)
                    hold on
                    x=(k-SPAN/2+SPAN*(c-1)/45+SPAN*mod(L-1,2)/90)*[1 1];
                    y=[starts_(SP) ends_(SP)];
                    if abs(y(2)-y(1))>4
                        continue
                    end
                    plot(x,y,'color',cols(L,:),'linewidth',.05)

                    x=(-SPAN/2+SPAN*(45*(k-1)+c-1)/(45*8)+SPAN*mod(L-1,2)/(45*8*2))*[1 1];
                    plot(x,y,'color',cols(L,:),'linewidth',.05)
                end
            end
        end
    end
    nlumps(:,k)=[trials.nlumps]';
    groups(:,k)=ones(45,1)*k;
end


subplot(5,3,1:3)
set(gca,'xtick',[])
ylabel('Fractional Cumulative Distance')

subplot(5,3,4:6)
set(gca,'xtick',[])
ylim([0 1.5])
ylabel('Fractional Distance to Target')

subplot(5,3,7:9)
set(gca,'xtick',[0:8])
set(gca,'xticklabel',{'All','1','2','3','4','5','6','7','8'});

ylabel('Time, sec')

set(gcf,'Color',backgray*[1 1 1])

subplot(5,3,14)
groups2=[groups(:); 0*groups(:)];
nlumps2=[nlumps(:); nlumps(:)];
boxplot(nlumps2,groups2,'orientation','vertical','notch','on','color','k','symbol','k+')
set(gca,'xtick',[1:9])
set(gca,'xticklabel',{'All','1','2','3','4','5','6','7','8'});
ylabel('Submovements Detected')
set(gca,'Color',backgray*[1 1 1])

labels={'Fractional Cumulative Distance','Fractional Distance to Target','Time, sec'};
group=[groupCell{:}];
group=[group group+5];
star=[starts{:}];
en=[ends{:}];

colorme=zeros(10,3);
for k=1:5
    colorme(k,:)=cols(k,:);
    colorme(k+5,:)=cols(k,:);
end

poses=zeros(10,1);
for k=1:5
    poses(k)=k-.11;
    poses(k+5)=k+.11;
end
    

for MET=1:3
    subplot(5,3,9+MET)
    cla
    dat=[star(MET,:), en(MET,:)];
    dat=real(dat);
    f=find((dat<3)&(dat>=0));
    dat=dat(f);
    gdat=group(f);
    boxplot(dat,gdat,'orientation','vertical','notch','on','color',colorme,'positions',poses,'width',.2,'symbol','k+')
    ylabel(labels{MET})
    set(gca,'xtick',1:5)
    set(gca,'xticklabel',{'1','2','3','4','5'})
    set(gca,'Color',backgray*[1 1 1])
end


