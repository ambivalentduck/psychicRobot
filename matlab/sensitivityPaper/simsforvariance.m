function simsforvariance(name,t,xvaf)

%% Step 1: Set up Monte Carlo Variance estimation

dat=paramsPopulator('burdet');
f=find(dat(:,4));
lf=length(f);

p=sobolset(2*lf,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

varied=p(1:1000,:); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.
%The number after p controls the number of points. Remember, this N*20 is the number of individual sims done.

%Clamp to +/- 3 stdevs to avoid support issues like negative segment masses
varied(varied>.999)=.999;
varied(varied<.001)=.001;

A=varied(:,1:lf);
B=varied(:,lf+1:2*lf);

AB=zeros(size(A,1),size(A,2),lf);
for k=1:lf
    AB(:,:,k)=A;
    AB(:,k,k)=B(:,k);
end

%% Step 2: Actually do the estimation

simA=repeatedSim(paramsPopulator(A),t,xvaf,'reflex',0,lf+1);
simB=repeatedSim(paramsPopulator(B),t,xvaf,'reflex',0,lf);

for k=1:lf
    k
    if k==1
        simAB=repeatedSim(paramsPopulator(AB(:,:,k)),t,xvaf,'reflex',0,lf-k);
    else
        simAB(k,:)=repeatedSim(paramsPopulator(AB(:,:,k)),t,xvaf,'reflex',0,lf-k);
    end
end
save(['BATCH',name,'.mat'],'simA','simB','simAB','t','xvaf')