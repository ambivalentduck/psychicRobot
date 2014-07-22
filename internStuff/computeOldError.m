clc
clear all

% Ss=[2 3 4 5 6 7 8 9];
Ss=2:9;
% Ss=23:30;

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
        if k >= 3*96+1 && k <= 4*96
%               In the intent phase(phase 3) the intent is also the cursor 
            [hand(k,S,1), hand(k,S,2), hand(k,S,3)]=maxperpendicular(trials(k).x,x0,x1);
            [cursor(k,S,1), cursor(k,S,2), cursor(k,S,3)]=maxperpendicular(trials(k).y,x0,x1);
            [intent(k,S,1), intent(k,S,2), intent(k,S,3)]=maxperpendicular(trials(k).y,x0,x1);
            fprintf('gg!!gg!!gg!!gg!!gg!!gg!!gg!!gg!!gg! %f !gg!!gg!!gg!!gg!!gg!!gg!!gg!!gg!!gg!!\n',k)
        else
%               In all the other phases (1,2,4, and 5) the cursor is the
%               hand
            [hand(k,S,1), hand(k,S,2), hand(k,S,3)]=maxperpendicular(trials(k).x,x0,x1);
            [cursor(k,S,1), cursor(k,S,2), cursor(k,S,3)]=maxperpendicular(trials(k).x,x0,x1);
            [intent(k,S,1), intent(k,S,2), intent(k,S,3)]=maxperpendicular(trials(k).y,x0,x1);
            fprintf(':D :D :D :D :D :D :D :D :D :D :D %f :D :D :D :D :D :D :D :D :D :D :D :D :D\n', k) 
        end
        [force(k,S,1), force(k,S,2), force(k,S,3)]=maxperpendicular(trials(k).f,x0,x1);

        onset=find(vecmag(trials(k).v)>.05,1,'first');
        term=find(vecmag(trials(k).cursor-ones(length(trials(k).t),1)*trials(k).targ)<.02,1,'first');
        reachT(k,S)=trials(k).t(term)-trials(k).t(onset);
        fprintf('lol %f lol\n', k)
    end
    toc
end

save('./Data/olderrors.mat','hand','intent','cursor','force','reachT')