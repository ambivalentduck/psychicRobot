clc
clear all

M=1000;

minN=2; %We see as few as two submotions 
meanN=4; %But 4 is the average
mu=meanN-minN;
dCmin=0; %.05^2; % Empirically, submotions are separated by about 50 ms
lambda=.15^2; % The average submotion is 200 ms

n=minN+floor(exprnd(mu,M,1));

D=zeros(M,1);
for k=1:M
    dCs2=exprnd(lambda,1,n(k));

    D(k)=sum(sqrt(dCmin+dCs2));
end

figure(5)
clf
hold on
[f,x]=ecdf(D);
plot(x,f,'b.')

%hist(Sn2)
[shift,nf,T]=fitShiftedGam(D)
plot(x,gamcdf(x-shift,nf,T),'r')
a=unique(n);
for k=a'
    i=find(n==k);
    [f,x]=ecdf(D(i));
    plot(x,f/(M/length(i)),'k')
end

xlim([0 1.5])
xlabel('Movement Duration, seconds')
