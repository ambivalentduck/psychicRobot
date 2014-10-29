clc
clear all

clin(2).bb=[38 39];
clin(2).fm=[43 45];
clin(2).wmf=[3.25 3.35];

clin(7).bb=[40 41];
clin(7).fm=[37 36];
clin(7).wmf=[3.14 2.88];

clin(8).bb=[35 37];
clin(8).fm=[40 38];
clin(8).wmf=[2.76 2.43];

subs=[2 7 8];
days=[1 6];
ks=[1 3];

for S=1:length(subs)
    stringS=num2str(subs(S),'%1.2i');
    for D=1:2
        day=days(D);
        %for k=1:3
        k=ks(D);

        load(['./org_data/outputsS',stringS,'D1F',num2str(k),'_org_session.mat']);

        p=Smoothed_pos;
        v=.01*[gradient(p(:,1)) gradient(p(:,2)) gradient(p(:,3))];
        a=.01*[gradient(v(:,1)) gradient(v(:,2)) gradient(v(:,3))];
        P=abs(dot(v',a'));
        [counts,centers]=hist(P,30);
        counts=counts/sum(counts);
        nzcounts=counts(counts~=0)';
        lognzcounts=log(nzcounts);
        nzcenters=centers(counts~=0)';
        Wz=[lognzcounts ones(size(lognzcounts))]\-nzcenters;
        R=cov(lognzcounts,nzcenters)/(std(lognzcounts)*std(nzcenters));
        R2=R(1,2)^2;
        
        p=vertcat(Pointer_arm.smooth_pos);
        v=.01*[gradient(p(:,1)) gradient(p(:,2)) gradient(p(:,3))];
        a=.01*[gradient(v(:,1)) gradient(v(:,2)) gradient(v(:,3))];
        P=abs(dot(v',a'));
        [counts,centers]=hist(P,30);
        nzcounts=counts(counts~=0)';
        lognzcounts=log(nzcounts);
        nzcenters=centers(counts~=0)';
        uWz=[lognzcounts ones(size(lognzcounts))]\-nzcenters;
        uR=cov(lognzcounts,nzcenters)/(std(lognzcounts)*std(nzcenters));
        uR2=uR(1,2)^2;
        
        data(S,D,k).W=Wz(1);
        data(S,D,k).R2=R2;
        data(S,D,k).uW=uWz(1);
        data(S,D,k).uR2=uR2;
        data(S,D,k).bb=clin(subs(S)).bb(D);
        data(S,D,k).fm=clin(subs(S)).fm(D);
        data(S,D,k).wmf=clin(subs(S)).wmf(D);
    end
end


W=[data.W];
R2=[data.R2];
uW=[data.uW];
uR2=[data.uR2];
bb=[data.bb];
fm=[data.fm];
wmf=[data.wmf];

movemeas=[R2' W'];
umovemeas=[uR2' uW'];

movelabs={'Percent Explained by Boltzmann','Wattage'};
clinmeas=[bb' fm' wmf'];
clinlabs={'Box+Blocks','Fugl-Meyer','Wolf Motor Function'};

figure(1)
clf
for k=1:2
    for kk=1:3
        subplot(2,3,3*(k-1)+kk)
        hold on
        plot(movemeas(:,k),clinmeas(:,kk),'r.')
        plot(umovemeas(:,k),clinmeas(:,kk),'k.')
        xlabel(movelabs{k})
        ylabel(clinlabs{kk})
    end
end

save('forYaz.mat','uW','uR2','W','R2','bb')