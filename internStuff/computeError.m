clc
clear all

Ss=23:30;

hand=zeros(5*96,8,3);
intent=zeros(5*96,8,3);
cursor=zeros(5*96,8,3);
force=zeros(5*96,8,3);
reachT=zeros(5*96,8);


tic
for S=1:8
    load(['./Data/output',num2str(Ss(S)),'.mat'])

    for k=1:5*96
        x0=trials(k).orig;
        x1=trials(k).targ;
        [hand(k,S,1), hand(k,S,2), hand(k,S,3)]=maxperpendicular(trials(k).x,x0,x1);
        [cursor(k,S,1), cursor(k,S,2), cursor(k,S,3)]=maxperpendicular(trials(k).cursor,x0,x1);
        [intent(k,S,1), intent(k,S,2), intent(k,S,3)]=maxperpendicular(trials(k).y,x0,x1);
        [force(k,S,1), force(k,S,2), force(k,S,3)]=maxperpendicular(trials(k).f,x0,x1);
        reachT(k,S)=trials(k).t(end)-trials(k).t(1);
    end
    toc
end

save('./Data/errors.mat','hand','intent','cursor','force','reachT')