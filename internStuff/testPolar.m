clc
clear all

nsubs=8;
exp1=2:9;
exp2=23:30;
subs=[exp1 exp2];
tlist=(1:480)';
b=floor((tlist-1)/96)+1;

%It's possible to start 2cm away from R=0
Rs=logspace(log10(2),log10(15),20);
thetabins=linspace(0,2*pi,16*3); %3 targs, 16 reaches per targ.

figure(1)
clf

colors=.2+.6*rand(8,3);

largest=0;

for S=1:8
    if S<20
        prefix='intern';
    else
        prefix='output';
    end
    load(['./Data/',prefix,num2str(S),'.mat'])
    ind=find([trials.targcat]==0,1,'first');
    center=trials(ind).targ;
    for targs=1:3
        indt=find([trials.targcat]==targs,1,'first');
        [theta,r]=cart2pol(trials(indt).targ(1)-center(1),trials(indt).targ(2)-center(2));
        [targs theta*180/pi+180+30]
        %With this transformation, theta*180/pi+180+30, you should be able
        %to SORT by transTheta and get an ordered list of targets. BUT it
        %may make more sense to replace that 30 with the most extreme 1 as
        %the wrap point.
    end
    
    thetas=zeros(16*3*4+16*3-1,length(Rs));
    thetasy=zeros(16*3*4+16*3-1,length(Rs));
    for k=1:length(trials)
        if trials(k).targcat==0
            continue
        end
        
        %Make a trials by distances matrix of thetas for hand.
        [trials(k).theta,trials(k).r]=cart2pol(trials(k).x(:,1)-center(1),trials(k).x(:,2)-center(2));
        %[trials(k).thetaX,trials(k).r]=cart2pol(trials(k).x(:,1)-center(1),trials(k).x(:,2)-center(2));
        %[trials(k).theta,trials(k).rV]=cart2pol(trials(k).v(:,1),trials(k).v(:,2));
        trials(k).r=trials(k).r*100;
        trials(k).theta=trials(k).theta*180/pi+180+30;
        for r_iter=1:length(Rs)
            [val,ind]=min(abs(trials(k).r-Rs(r_iter)));
            thetas(k,r_iter)=trials(k).theta(ind);
        end
        
        %Make a trials by distances matrix of thetas for intent.
        [trials(k).thetay,trials(k).ry]=cart2pol(trials(k).y(:,1)-center(1),trials(k).y(:,2)-center(2));
        %[trials(k).thetayX,trials(k).ry]=cart2pol(trials(k).y(:,1)-center(1),trials(k).y(:,2)-center(2));
        %[trials(k).thetay,trials(k).ryV]=cart2pol(gradient(trials(k).y(:,1)),gradient(trials(k).y(:,2)));
        trials(k).ry=trials(k).ry*100;
        trials(k).thetay=trials(k).thetay*180/pi+180+30;
        for r_iter=1:length(Rs)
            [val,ind]=min(abs(trials(k).ry-Rs(r_iter)));
            thetasy(k,r_iter)=trials(k).thetay(ind);
        end
        
    end
    
    targs=[trials.targcat]';
    dummytargs=1+(floor((0:47)/16))';
    boolNotZero=targs~=0;
    
    Hx=log2(48)-(1/48)*3*16*log2(16);
    Hy=Hx;
    
    effiencies=zeros(6,length(Rs));
    
    for block=1:6
        if block~=6
            f=find((b==block)&boolNotZero);
        else
            f=find((b==3)&boolNotZero);
        end
        targlabs=targs(f);
        for r_iter=1:length(Rs)

            if block~=6
                [sorted,inds]=sort(thetas(f,r_iter));
            else
                [sorted,inds]=sort(thetasy(f,r_iter));
            end
            sortedTarglabs=targlabs(inds);

            Hxy=log2(48);
            for c=1:3
                for cc=1:3
                    Nij=sum((sortedTarglabs==c)&(dummytargs==cc));
                    if Nij~=0
                        Hxy=Hxy-(1/48)*Nij*log2(Nij);
                    end
                end
            end
            efficiencies(block,r_iter)=(Hx+Hy-Hxy)/Hx;
        end
    end
    efficiencies'
    fade=.8;
    handfade=[fade fade 1];
    intentfade=[1 fade fade];
    subplot(3,1,1)
    hold on
    plot(Rs,efficiencies(2,:),'color',handfade)
    subplot(3,1,2)
    hold on
    plot(Rs,efficiencies(3,:),'color',handfade)
    plot(Rs,efficiencies(6,:),'color',intentfade)
    subplot(3,1,3)
    hold on
    plot(Rs,efficiencies(4,:),'color',handfade)
    
    xsubEff(:,:,S)=efficiencies;
end

figure(1)
hand=[0 0 .5];
subplot(3,1,1)
plot(Rs,mean(xsubEff(2,:,:),3),'color',hand,'linewidth',2)
xlim([2 15])
ylim([0 1.05])
title('Block 2')


subplot(3,1,2)
plot(Rs,mean(xsubEff(3,:,:),3),'color',hand,'linewidth',2)
plot(Rs,mean(xsubEff(6,:,:),3),'color',[.5 0 0],'linewidth',2)
xlim([2 15])
ylim([0 1.05])
ylabel('Visual-Motor Efficiency')
title('Block 3')

subplot(3,1,3)
plot(Rs,mean(xsubEff(4,:,:),3),'color',hand,'linewidth',2)
xlim([2 15])
ylim([0 1.05])
xlabel('Progress, cm')
title('Block 4')


%% Make fig 2
figure(2)
clf
hold on
alpha=.15;
m2=mean(xsubEff(2,:,:),3);
s2=std(xsubEff(2,:,:),0,3);
m3=mean(xsubEff(3,:,:),3);
s3=std(xsubEff(3,:,:),0,3);
m4=mean(xsubEff(4,:,:),3);
s4=std(xsubEff(4,:,:),0,3);
m6=mean(xsubEff(6,:,:),3);
s6=std(xsubEff(6,:,:),0,3);

Rs=Rs-2;
c=1;
while ttest(squeeze(xsubEff(6,c,:)-xsubEff(2,c,:)))
    c=c+1;
end
Rs(c)

h=fill([Rs wrev(Rs)],[m2+s2 wrev(m2-s2)],'g');
set(h,'edgealpha',0,'facealpha',alpha);
h2=plot(Rs,m2,'color','g','linewidth',2);

h=fill([Rs wrev(Rs)],[m3+s3 wrev(m3-s3)],[0 0 .7]);
set(h,'edgealpha',0,'facealpha',alpha);
h3=plot(Rs,m3,'color',[0 0 .7],'linewidth',2);

h=fill([Rs wrev(Rs)],[m4+s4 wrev(m4-s4)],'c');
set(h,'edgealpha',0,'facealpha',alpha);
h4=plot(Rs,m4,'color','c','linewidth',2);

h=fill([Rs wrev(Rs)],[m6+s6 wrev(m6-s6)],[.7 0 0]);
set(h,'edgealpha',0,'facealpha',alpha);
h6=plot(Rs,m6,'color',[.7 0 0],'linewidth',2);

xlim([0 13])
ylim([0 1.05])
legend([h2 h3 h6 h4],{'Block 2, Hand','Block 3, Hand','Block 3, Intent','Block 4, Hand'})
xlabel('Progress, cm')
ylabel('Visual-Motor Efficiency')
title('Across Subject Averages')