clc
clear all

pams=zeros(8,2);

for k=1:8
    load(['../Data/curlkick/curlkick',num2str(k),'.mat'])
    [pams(k,:),dmse(k)]=correctcdf(trials,k)
end

inds=(1:8)~=3;

signtest(dmse(inds))
mu=mean(pams(inds,:))
mu(1)
sd=std(pams(inds,:))
sd(1)
    