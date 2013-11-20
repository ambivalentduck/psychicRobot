clc
clear all

submax=8;
backgray=.6;

starts={};
ends={};

figure(666)
clf
subplot(5,3,14)
nlumps=zeros(45,submax);
groups=nlumps;
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
        starts{k,c}=zeros(length(f2),3);
        ends{k,c}=zeros(length(f2),3);
        for cc=1:length(f2)
            kk=f(f2(cc));
            inds=trials(kk).lumps(c).inds;
            rat=targdists(targetcat(kk));
            starts{k,c}(cc,1:3)=[trials(kk).cumdist(inds(1))/rat trials(kk).yurot(inds(1),1)/rat trials(kk).tcat(inds(1))];
            ends{k,c}(cc,1:3)=[trials(kk).cumdist(inds(end))/rat trials(kk).yurot(inds(end),1)/rat trials(kk).tcat(inds(end))];
        end
    end

    nlumps(:,k)=[trials.nlumps]';
    groups(:,k)=ones(45,1)*k;
    names{k+1}=num2str(k);
end

subnums={};
for k=1:8
    for c=1:5
        starts{k,c}=starts{k,c}';
        ends{k,c}=ends{k,c}';
        subnums{k,c}=k*ones(size(starts{k,c}));
    end
end

names{1}='All';
groups2=[groups(:); 0*groups(:)];
nlumps2=[nlumps(:); nlumps(:)];
%boxplot(nlumps2,groups2,'orientation','horizontal','notch','on','labels',names)
boxplot(nlumps2,groups2,'orientation','vertical','notch','on','labels',names)
%boxplot(nlumps2,groups2,'orientation','vertical','plotstyle','compact','labels',names)
xlabel('Subject Number')
ylabel('Submovements Detected')
set(gca,'Color',backgray*[1 1 1])

yl=[3 1.5 2.5];
tits={'Cumulative Distance Ratio','Progress Ratio','Time, s'};
pos=0:8;
pos=[pos-.1 pos+.1];

cols=zeros(18,3);
% for k=1:18
%     if mod(k,2)
%         cols(k,:)=[0 0 1];
%     else
%         cols(k,:)=[1 0 0];
%     end
% end
cols(1:9,3)=1;
cols(10:18,1)=1;

for k=1:4
    dats=[starts{:,k}]';
    date=[ends{:,k}]';
    sms=[subnums{:,k}]';
    for met=1:3
        subplot(5,3,met+3*(k-1))
        dat2=[dats(:,met);dats(:,met);date(:,met);date(:,met)];
        f=find((~imag(dat2))&(dat2<50));
        sms2=[sms(:,met);0*sms(:,met);sms(:,met)+9;0*sms(:,met)+9];
        boxplot(dat2(f),sms2(f),'orientation','vertical','plotstyle','compact','positions',pos,'colors',cols)
        if k==1
            title(tits{met})
        end
        if met==1
            ylabel(['Submovement ',num2str(k)])
        else
            ylabel('')
        end
        ylim([0 yl(met)])
        set(gca,'xtick',0:8)
        set(gca,'xticklabel',names)
        set(gca,'Color',backgray*[1 1 1])
    end
end

set(gcf,'Color',backgray*[1 1 1])




