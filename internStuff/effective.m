clc
clear all

%cleanup

global x0 l1 l2

subs=[2:9 23:30];
center=[0 .48];
kp1=[15, 6;6, 16];
kp2=[8, 2;2, 5];

P=0:5:360;
eps=1/1000;

modelstiff=zeros(max(subs),length(P));
Smodelstiff=zeros(max(subs),length(P));

for S=subs
    if S<20
        prefix='intern';
    else
        prefix='output';
    end
    load(['./Data/',prefix,num2str(S),'.mat'])
    x0=params.shoulder;
    l1=params.l1;
    l2=params.l2;
    hash=num2str(S);
    aName=['getAlpha',hash];
    fName=['fJ',hash];
    if ~exist([fName,'.m'],'file')
        makeJacobians(aName,fName);
        continue
    end
    fJ=str2func(fName);
    
    q=ikin(center);
    for p=1:length(P)
        qp=ikin(center+eps*[sind(P(p)) cosd(P(p))]);
        if S<15
            tau=kp1*(q-qp);
            Jt=fJ(qp)';
            modelstiff(S,p)=norm(Jt\tau)/eps;
            Smodelstiff(S,p)=norm(Jt\tau)/eps;
        else
            tau=kp1*(q-qp);
            Jt=fJ(qp)';
            Smodelstiff(S,p)=norm(Jt\tau)/eps;
            
            tau=kp2*(q-qp);
            Jt=fJ(qp)';
            modelstiff(S,p)=norm(Jt\tau)/eps;
        end
        
    end
end

figure(299)
clf
hold on
for S=subs
    plot(S*ones(length(P),1)-.1,modelstiff(S,:),'b.')
end

save('modelstiff.mat','modelstiff','Smodelstiff')
