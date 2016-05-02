clc
clear all

%% Get empirical cdf
for SUBNUM=7 %1:8
    resultsMat(1,SUBNUM)=SUBNUM;
    
    load(['../Data/curlkick/curlkick',num2str(SUBNUM),'.mat'])
    
    tcats=[trials.targetcat];
    dcats=[trials.disturbcat];
    
    f=find((tcats~=3)&~dcats);
    Vmax=zeros(size(f));
    
    for k=1:length(f)
        tr=trials(f(k));
        Vmax(k)=max(vecmag(trials(f(k)).v)).^2;
    end
    
    figure(SUBNUM)
    clf
    hold on
    ecdf(Vmax,'bounds','on')
    mu=mean(Vmax);
    sigma=std(Vmax);
    [f,x]=ecdf(Vmax);
    hgauss=plot(x,normcdf(x,mu,sigma),'r')
    resultsMat(2,SUBNUM)=jbtest(Vmax)
    
    %% Max n-based cdf
    
    nMin=2;
    thetas=[1.93, 1.8, 1.55,1.18,1.4,1.49,2.1,1.42];
    theta=thetas(SUBNUM);
    M=10000;
    
    ns=nMin+round(exprnd(theta,M,1));
    
    % [f,x]=ecdf(ns);
    % figure(1)
    % clf
    % plot(x(1:end-1),log(1-f(1:end-1)))
    % p=1./polyfit(x(1:end-1),log(1-f(1:end-1)),1)
    
    maxes=zeros(M,1);
    for k=1:M
        %vals(k).x=rand(1,ns(k));
        vals(k).x=exprnd(1,1,ns(k));
        vals(k).x=vals(k).x/sum(vals(k).x);
        maxes(k)=max(vals(k).x);
    end
    
    %maxexp=@(params) (1-exp(-(x-params(3)).^2/params(1))).^params(2);
    %maxexp=@(params) (1-exp(-x.^2/params(1))).^params(2);
    [fm,xm]=ecdf(maxes);
    
    relate=@(params) interp1(params(1)*xm(2:end)+params(2),fm(2:end),x,'nearest','extrap');
    
    minme=@(params) sum((f-relate(params)).^2);
    %minme=@(params) sum((f-maxexp(params)).^2);
    
    guess=[max(xm) 1; min(xm) 1]\[max(x);min(x)];
    %guess=[.1; .4];
    
    ps=fminsearch(minme,guess);
    resultsMat(3:4,SUBNUM)=ps;
    %ps=fminunc(minme,guess,optimoptions('fminunc','TolX',1e-12,'TolFun',1e-29,'MaxFunEvals',10000))
    mseFit=minme(ps)/length(x)
    mseGauss=mean((f-normcdf(x,mu,sigma)).^2)
    
    resultsMat(5:6,SUBNUM)=[mseFit; mseGauss];
    
    hfit=plot(x,relate(ps),'g')
end

resultsMat
[win,howmuch]=ttest(resultsMat(5,:)-resultsMat(6,:))

%% Polish figure
hfig=figure(7)
ylabel('Cumulative Probability')
xlabel('Peak Speed, m/s')

hcdf=findall(hfig,'Type','line')

legend(hcdf([5 2 1]),{'Data','Gaussian','Model'})
%xlim([.1 .3])
