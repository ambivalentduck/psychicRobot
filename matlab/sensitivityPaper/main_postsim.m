%Just in case...
save('all_kick_workspace.mat')


%% Burdet first


bins=0:.005:.15;

[trash,ref]=getMUE(bins,0*bins,yex);
numerical_accuracy_burdet=getMUE(bins,ref,yex)

for k=1:size(OAT,1)
    for kk=1:size(OAT,2)
        OAT(k,kk).mue=getMUE(bins,ref,OAT(k,kk).y);
    end
end

mue=[OAT.mue];
mue=reshape(mue,size(OAT,1),size(OAT,2))

numP=size(simAB,1);
N=length(simA);
for k=1:N
    simA(k).mue=getMUE(bins,ref,simA(k).y);
    simB(k).mue=getMUE(bins,ref,simB(k).y);
    for kk=1:numP
        simAB(kk,k).mue=getMUE(bins,ref,simAB(kk,k).y);
    end
end

mueA=[simA.mue];
mueB=[simB.mue];
mueAB=[simAB.mue];
mueAB=reshape(mueAB,size(simAB,1),size(simAB,2));

varY=var(mueAB(:))

S=zeros(numP,1);
ST=S;

for k=1:numP
    S(k)=1/(2*N)*sum((mueA-mueAB(k,:)).^2);
    ST(k)=1/N*sum(mueB.*(mueAB(k,:)-mueA));
end
S=S/varY
ST=ST/varY

dat=paramsPopulator('burdet');
figure(6)
clf
barh(1:numP,[S ST])
set(gca,'ytick',1:numP)
set(gca,'yticklabel',names(dat(:,4)>0))

%% Then Shadmuss

[trash,refsm]=getMUE(bins,0*bins,yexsm);
numerical_accuracy_shadmuss=getMUE(bins,refsm,yexsm)

for k=1:size(OATSM,1)
    for kk=1:size(OATSM,2)
        OATSM(k,kk).mue=getMUE(bins,refsm,OATSM(k,kk).y);
    end
end

muesm=[OATSM.mue];
muesm=reshape(muesm,size(OATSM,1),size(OATSM,2))




