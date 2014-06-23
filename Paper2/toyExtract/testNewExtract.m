clc
clear all

global M B measuredVals measuredTime

%% Constant K, pulse Force
tspan=.5;
tstep=.005;
ta=(0:tstep:tspan)/tspan;

M=.41;
B=2.3;

xi=[.1; 0];
xf=[0; 0];

xd=(xi*ones(size(ta))+(xf-xi)*(10*ta.^3-15*ta.^4+6*ta.^5))';
vd=((xf-xi)*(30*ta.^2-60*ta.^3+30*ta.^4)/tspan)';
ad=((xf-xi)*(60*ta-180*ta.^2+120*ta.^3)/(tspan^2))';
% tshort=0:.001:tspan;
% xd=([1; 0]*sin(tshort))';
% vd=([1; 0]*cos(tshort))';
% ad=([1; 0]*-sin(tshort))';

t=0:tstep:1.5;
Y=[[xd vd ad];zeros(length(t)-length(ta),6)];
Y(:,1)=Y(:,1)+1;
F=zeros(length(t),2);
F((t>.1)&(t<.2),1)=-10;
%F(:,1)=5*cos(3*t);
K=16*ones(length(t),2);

Xh=forward(t,Y,F,K);

%For shits and giggles, make a bank
nBank=8;
maxRatio=.4;
minK=16;
bank=zeros(length(t),nBank);
eps=(1+maxRatio*((1:nBank)-1)/(nBank-1))*minK
for k=1:nBank
    Xtemp=reverse(t,Xh,F,eps(k)*ones(length(t),2));
    bank(:,k)=Xtemp(:,1);
    if k==1
        %Want the lowest suspected overestimate in place
        [T,Fcorr]=ode45(@integrateForce,t,zeros(4,1));
    end
end

threshFcorr=.05*max(abs(Fcorr(:,1))); % 5% of a 1N pulse is nothing.

minned=zeros(length(t),2);
for k=1:length(t)
    warning off
    minme=@(m) std(bank(k,:)-m.*eps.*Fcorr(k,1));
    minned(k,1)=fminunc(minme,1,optimset('display','off'));
end

inds=find(abs(Fcorr(:,1))>threshFcorr);
andmetoo=@(m) sum(corr(Fcorr(inds,1),bank(inds,1)-minned(inds,1).*(eps(1)-m).*Fcorr(inds,1)).^2);
%andmetoo=@(m) sum(corr(Fcorr(inds,1),bank(inds,1)-minned(inds,1).*(eps(1)-m).*Fcorr(inds,1)).^2);
minned(:,2)=fminunc(andmetoo,5,optimset('display','off'));

figure(1)
clf
hold on
for k=1:nBank
    plot(t,bank(:,k)-minned(:,1).*(eps(k)-minned(:,2)).*Fcorr(:,1),'r')
    plot(t,bank(:,k)-minned(:,1).*(eps(k)-16).*Fcorr(:,1),'r')
end
plot(t,10*Fcorr(:,1)+1,'b')
plot(t,Y(:,1),'g')
title('Constant K by minimizing Pearson''s R^2')
minned(1,2)

figure(2)
clf
hold on
mY=norm(Y(:,1));
corrFcorr=minned(:,1).*Fcorr(:,1);
corrFcorr=corrFcorr-mean(corrFcorr); %have to make stationary
mFcorr=dot(corrFcorr,corrFcorr); %norm(corrFcorr);
for k=1:nBank
    tempcalc=bank(:,k)-eps(k)*minned(:,1).*Fcorr(:,1);
    tempcalc=tempcalc-mean(tempcalc);
    d=dot(tempcalc,corrFcorr)/(mFcorr);
    ycalc=bank(:,k)-(eps(k)+d)*corrFcorr;
    plot(t,ycalc-mean(ycalc),'r')
    maxerror=max(abs(ycalc-Y(:,1)))
end
plot(t,Y(:,1)-mean(Y(:,1)),'g')
plot(t,corrFcorr*10,'color',[.5 .5 .5])
title(['Local non-stationarity prevents global use of dot product. Kcalc=',num2str(-d)])
return

figure(3)
clf
hold on
corrFcorr=minned(:,1).*Fcorr(:,1);
Fspan=60;
kcalc=zeros(length(t),1);
kcalc(1:Fspan)=minK;
for k=1:nBank
    %Now have to have yet another for-loop to deal with subindexing.
    %Locally stationary we can probably do, so this is in many/most senses
    %ideal for time-varying K.
    for kk=Fspan+1:length(t)
        inds=max(1,kk-Fspan):kk;

        ibank=bank(inds,k);

        iFc=corrFcorr(inds);
        iFc=iFc-mean(iFc);
        mFcorr=norm(iFc);

        d2=inf;
        dcum=0;
        its=0;
        while (abs(d2)>1e-3)&&(its<15)
            its=its+1;
            tempcalc=ibank-(eps(k)+dcum).*iFc;
            tempcalc=tempcalc-mean(tempcalc);
            d=dot(tempcalc,iFc)/(mFcorr); %*norm(bank(:,k)-eps(k)*corrFcorr))
            dcum=dcum+d;
            if its==1
                dcum
            end
            ycalc=ibank-(eps(k)+dcum).*iFc;
            ycalc=ycalc-mean(ycalc);
            d2=dot(iFc,ycalc)/(mFcorr*norm(ycalc));
        end
        kcalc(kk)=-dcum;
        [k kk]
    end
    plot(t,bank(:,k)-(eps(k)-kcalc).*corrFcorr,'r.')
    plot(t,kcalc/10,'y')
    plot(t,bank(:,k)-(eps(k)-0)*corrFcorr,'c')
    drawnow
end
plot(t,Y(:,1),'g')
plot(t,corrFcorr*100,'color',[.5 .5 .5])
title('We''d like to make everything local anyways, but does it work?')


% figure(2)
% clf
% hold on
% for G=-15:3:15
% for k=1:nBank
%     plot(t,gradient(bank(:,k)-minned(:,1).*(eps(k)-minned(:,2)+G).*Fcorr(:,1)),'r')
%     %plot(t,bank(:,k)-minned(:,1).*(eps(k)-minned(:,2)).*Fcorr(:,1),'r')
% end
% end
% plot(t,gradient(10*Fcorr(:,1)+1),'b')
% plot(t,gradient(Y(:,1)),'g')
% title('Constant K')


return

%% Time-varying K, pulse Force

filtN=10;
K=15+5*filtfilt(ones(filtN,1)/filtN,1,randn(length(t)+2*filtN,1));
K=K(1+filtN:end-filtN); %Has quirks starting and stopping
K=[K K];

figure(2)
clf
hold on
plot(t,K(:,1)/10,'color',[.5 .5 .5])

Xh=forward(t,Y,F,K);

plot(t,Xh(:,1),'k')

%For shits and giggles, make a bank
nBank=8;
maxRatio=.4;
minK=15.5;
bank=zeros(length(t),nBank);
eps=(1+maxRatio*((1:nBank)-1)/(nBank-1))*minK
for k=1:nBank
    Xtemp=reverse(t,Xh,F,eps(k)*ones(length(t),2));
    bank(:,k)=Xtemp(:,1);
end

K0=15;
minned=zeros(length(t),2);
Fspan=10; %50 ms isn't awful
%threshFcorr=.05*max(abs(Fcorr(:,1))); % 5% of a 1N pulse is nothing.
threshFcorr=6.1428e-4; %Hardcoded version of the above. A tiny number.

figure(3)
clf
hold on

%The bank can be constant, but Fcorr ideally is not and needn't be anyways.
for k=1:length(t)
    %Use the last known K to update Fcorr, which. It conveniently always starts at or near zero.
    if (k==1) %||1
        Fcorr(1,:)=[0 0 0 0];
        measuredVals(1,9:10)=K0;
    else
        measuredVals(k,9:10)=minned(k-1,2); %Have to use a one-sample delay.
        [T,FcorrSim]=ode45(@integrateForce,t(k-1:k),Fcorr(k-1,:));
        Fcorr(k,:)=FcorrSim(end,:);
    end
    %Calculate the local epsilon scaling factor
    minme=@(m) std(bank(k,:)-m.*eps.*Fcorr(k,1));
    minned(k,1)=fminunc(minme,1,optimset('display','off'));

    %Calculate a somewhat lagged version of K
    inds=max(1,k-Fspan):k;
    tempF=Fcorr(inds,1);
    subinds=find(abs(tempF)>threshFcorr);
    if k==1
        minned(k,2)=K0;
    elseif length(subinds)<5
        minned(k,2)=minned(k-1,2);
    else
        andmetoo=@(m) sum(corr(Fcorr(inds(subinds),1),bank(inds(subinds),:)-repmat(minned(inds(subinds),1),1,nBank).*(repmat(eps,length(subinds),1)-m).*repmat(Fcorr(inds(subinds),1),1,nBank)).^2);
        minned(k,2)=fminunc(andmetoo,5,optimset('display','off'));
        plot(t(inds(subinds)),Fcorr(inds(subinds),1),'.','color',[.5 .5 .5])
        plot(t(inds(subinds)),bank(inds(subinds),1)-minned(inds(subinds),1).*(eps(1)-minned(k,2)).*Fcorr(inds(subinds),1),'c.')
        drawnow
    end
end

figure(2)
f=find(abs(minned(:,2)-15)<=5);
plot(t(f),minned(f,2)/10,'.')

return

for k=1:nBank
    plot(t,bank(:,k)-minned(:,1).*(eps(k)-minned(:,2)).*Fcorr(:,1),'r')
    %plot(t,bank(:,k)-minned(:,1).*(eps(k)-minned(:,2)).*Fcorr(:,1),'r')
end
plot(t,10*Fcorr(:,1)+1,'b')
plot(t,Y(:,1),'g')

return

epsilon1=2;
epsilon2=2.1;
Ye1=reverse(t,Xh,F,epsilon1+K);

epsilon2=2.1;
Ye2=reverse(t,Xh,F,epsilon2+K);

figure(1)
clf
hold on
sc=2 %mean(abs(Ye1(:,1)-Fcorr(:,1)))*mean(abs(Fcorr(:,1)));
plot(t,Y(:,1),'g',t,Xh(:,1),'k',t,Ye1(:,1),'r',t,Ye1(:,1)-sc*Fcorr(:,1),'b')

return
inds=find(Fcorr(:,1)~=0);
%inds=inds(10:end);
figure(2)
clf
hold on
plot(t,Fcorr(:,1),'color',[.5 .5 .5])
plot(t,Ye2(:,1)-Ye1(:,1),'b')

% rSpan=4;
% for k=1:length(t)
%     inds=max(1,k-rSpan):min(length(t),k+rSpan);
%     m(k)=[Fcorr(inds,1)]\(Ye2(inds,1)-Ye1(inds,1));
% end




