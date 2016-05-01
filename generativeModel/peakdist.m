clc
clear all

load('../Data/curlkick/curlkick1g.mat')

tcats=[trials.targetcat];
dcats=[trials.disturbcat];

figure(1)
clf
hold on

colors='rgb';
dirs=[1 2 4];

for DERP=1 %1:3
    f=find((tcats==dirs(DERP))&~dcats);
    Vmax=zeros(size(f));
    
    for k=1:length(f)
        tr=trials(f(k));
        Vmax(k)=max(vecmag(trials(f(k)).v));
    end
end

figure(1)
clf
hold on
ecdf(Vmax,'bounds','on')
mu=mean(Vmax);
sigma=std(Vmax);
[f,x]=ecdf(Vmax);
plot(x,normcdf(x,mu,sigma),'r')
jbtest(Vmax)
jbtest(randn(1000,1))

%maxexp=@(params) (1-exp(-(x-params(3)).^2/params(1))).^params(2);
maxexp=@(params) (1-exp(-x.^2/params(1))).^params(2);


minme=@(params) sum((f-maxexp(params)).^2);

ps=fminunc(minme,[.0301 24],optimoptions('fminunc','TolX',1e-12,'MaxFunEvals',1000))

plot(x,maxexp(ps),'g')