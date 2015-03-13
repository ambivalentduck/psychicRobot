clc
clear all

load('free_exp_05stroke.mat')

en=dot(v',v');

figure(75)
clf
inds=1:10:size(x,1);
plot3(x(inds,1),x(inds,2),en(inds),'.','markersize',1)

figure(76)
clf
upper=prctile(x,99);
lower=prctile(x,1);

binsx=[linspace(lower(1),upper(1),64)];
binsy=[linspace(lower(2),upper(2),64)];

counts=zeros(64,64);

for k=2:length(binsx)
    for kk=2:length(binsy)
        counts(k,kk)=sum((((x(:,1)>binsx(k-1))&(x(:,1)<=binsx(k))))&(((x(:,2)>binsy(kk-1))&(x(:,2)<=binsy(kk)))));
    end
end
       

%counts(counts==0)=1;
%counts=counts/sum(counts);

%surf(log(counts))
