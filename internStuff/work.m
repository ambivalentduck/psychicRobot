clc
clear all

load('./Data/errors.mat')
hand(:,:,2)=force(:,:,2)./hand(:,:,2);

Ss=23:30;
%Ss=1:8;

allT=1:5*96;
phases=zeros(5*96,6);
for k=[1 2 4 5]
    phases(96*(k-1)+1:k*96,k)=1;
end
phases=phases(:,[1 2 3 6 4 5]);
phases(2*96+1:2.5*96,3)=1;
phases(2.5*96+1:3*96,4)=1;

phases=phases==1;

workP=zeros(8*5*96,1);
stiffness=workP;
X=zeros(8*5*96,2);
for S=1:8
    load(['./Data/output',num2str(Ss(S)),'.mat'])
    
    inds=5*96*(S-1)+1:5*96*S;
    for k=1:length(trials)
        workP(inds(k))=sum(dot(trials(k).f,trials(k).v));
    end
    stiffness(inds)=hand(:,S,2);
    
    X(inds,1)=phases*(1:6)';
    X(inds,2)=S+inds*0;
end

figure(7)
[p,table,stats]=anovan(workP,X,'display','off')
multcompare(stats);

figure(8)
[p,table,stats]=anovan(workP./stiffness,X,'display','off')
multcompare(stats);
