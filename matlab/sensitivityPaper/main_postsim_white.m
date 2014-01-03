load all_white_workspace.mat

%% fig 1

qscale=.001;
SKIP=2;
figure(1)
clf
hold on
plot(x(:,1),x(:,2),'k.',xsim(:,1),xsim(:,2),'k-')
quiver(xsim(1:SKIP:end,1),xsim(1:SKIP:end,2),qscale*xvaf(1:SKIP:end,7),qscale*xvaf(1:SKIP:end,8),0,'Color',.6*[1 1 1])
plot(yex(:,1),yex(:,2),'ko')
plot(ycross(:,1),ycross(:,2),'k-.')
axis equal
axis off

bins=0:.005:.15;
[trash,ref]=getMUE(bins,0*bins,yex);
cross_accuracy=getMUE(bins,ref,ycross)

plot([0 .15],.486*[1 1],'k','linewidth',3)
text(0,.483,'15 cm')
plot(.0005*[1 1],[.505 .515],'Color',.6*[1 1 1],'linewidth',3)
text(-.001,.5075,'10 N','rotation',90)

ylim([0.47,0.53])
xlim([-.01,0.16])
return

%% fig 2

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
H=barh(1:numP,[S ST])
set(H(1),'facecolor',.2*[1 1 1])
set(H(2),'facecolor',.8*[1 1 1])
set(gca,'ytick',1:numP)
set(gca,'yticklabel',names(dat(:,4)>0))
xlabel('Fraction of Total Variance')
ylim([0 21])
xlim([-.05 1])

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




