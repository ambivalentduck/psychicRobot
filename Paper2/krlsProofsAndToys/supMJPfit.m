function [cost, grad]=supMJPfit(P)

global yfit x0 t dt

%% Step 0 P->w,ts

r=reshape(P,length(P)/3,3);
w=r(:,1:2);
ts=r(:,3);

regparam=1/1000;

%% Step 1 cost....5 sum error.^2 because why not? + Regularization

inds=ts>2*dt;

tinds=t(inds)';
ts2=ts(inds)/2;
uppers=tinds+ts2;
lowers=tinds-ts2;
[summed,kerns]=supMJP(x0,w(inds,:),lowers,uppers,t);
erroryfit=summed(:,3:4)-yfit;

overlap=zeros(length(ts),1);
gradw=zeros(length(ts),2);
gradts=zeros(length(ts),1);

kk=1:length(ts);
for k=kk(inds)
    overlap(k)=sum(max(0, min(uppers(k), uppers) - max(lowers(k), lowers))>0);
    gradw(k,:)=sum([erroryfit(:,1).*kerns(:,k) erroryfit(:,2).*kerns(:,k)]); %After summing, should be a row vector
    gradts(k)=sum(sum(erroryfit,2).*((45*(t' - t(k)).^2)/(ts(k)^4) - (150*(t' - t(k)).^4)/(ts(k)^6) - 15/(8*ts(k)^2))) + overlap(k)*regparam;
end
overlap=overlap-2; %Self and 1. Ie. Never more than 2-coverage.

cost=sum(sum(erroryfit.^2))+sum(overlap)*regparam;

grad=[gradw(:);gradts];

