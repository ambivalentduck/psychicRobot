clc
clear all

zed=0:.01:5;
lz=length(zed);
lambda=.0838; %Empirical, sub 1
theta=.01; %.5301; %Empirical, sub 1

maxn=6;
pon=exp(-(1:maxn)/theta);
%pon(1)=0;
pon=pon/sum(pon);

p=1/lambda*exp(-zed.^2/lambda)*2.*zed;
p=p/sum(p);
pgn=p;

pn=zeros(maxn,lz);

%% setup

%CDF = 1-exp(-y^2/lambda), no surprise there. Y = peak to peak duration
figure(5)
clf
subplot(1,2,1)
hold on

for k=1:maxn
    plot(zed,pgn(1:lz),'k')
    mu=sum(zed.*pgn(1:lz));
    sigma=sqrt(sum(zed.^2.*pgn(1:lz))-mu^2);
    npdf=normpdf(zed,mu,sigma);
    npdf=npdf/sum(npdf);
    plot(zed,npdf,'r')
    
    [mp,imp]=max(pgn);
    text(zed(imp),mp,num2str(k),'verticalalignment','bottom')
    
    f=find(pgn(1:lz)>.00001);
    cp=cumsum(pgn(f));
    rnd=rand(10000,1);
    samples=interp1(cp,zed(f),rnd);
    skewness(samples)
    pgn=conv(p,pgn);
    
    pn(k,:)=pgn(1:lz)'*pon(k);
end
title('Duration Depends on Submotion Count')
ylabel('Probability')
xlabel('Movement Duration')
legend('Model','Best-Fit Gaussian')

%% whole pdf

subplot(1,2,2)
hold on
spn=sum(pn);
plot(zed,cumsum(spn),'k')
%mu=sum(zed.*spn);
%sigma=sqrt(sum(zed.^2.*spn)-mu^2);
%ncdf=normcdf(zed,mu,sigma);
%plot(zed,ncdf,'r')

title('Summed Given Expected Counts')
xlabel('Movement Duration')
ylabel('Cumulative Probability')

%% add in real data

load('../Data/curlkick/curlkick1g.mat')

tcats=[trials.targetcat];
dcats=[trials.disturbcat];

f=find((tcats==4)&~dcats);

segments=[.4 .3 .2 .1];
S=zeros(length(f),length(segments));
RW=S;

for k=1:length(f)
    tr=trials(f(k));
    [S(k,:),RW(k,:)]=getTsMetric(tr.x,tr.v,tr.a,tr.t,tr.x(1,:),tr.x(end,:),segments);
end

ecdf(S(:,end),'bounds','on');

[mu,sigma]=normfit(S(:,end));
plot(zed,normcdf(zed,mu,sigma),'r')
xlim([.3 .6])
%% wtf

figure(6)
hist(S,50)

%% less than one

figure(7)
clf
hold on

fracsOfInt=[.9 1];

for each=fracsOfInt
    pInferred = ifft(fft(p).^each);
    
    plot(zed(1:10),real(pInferred(1:10)))
end