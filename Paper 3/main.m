clc
clear all

c=dir('./org_data');

%Specify subject numbers

%Each subject comes "with" a library of movements
%--Yazan owes Justin a means of getting this part right.
%--Reaches, to the center, with explicit encouragement to be bimanual with
%  NO external force disturbances

%For each movement, load its shit into a structure (and compute A), while doing that get Xrange

%Compute Xgrid

%For each movement, compute Ux, Ut

%From all movements, compile Ux as a statistical object

%Test P(x|movement)=P(x|U,T) for some T

%Use the confidence limit from the above in order to test Uhealthy=Uimpaired

%Plot voxels that fail the test as translucent cubes with spaghetti drawn on top
%--Maybe plot degree of difference as transparent-to-opaque


cMin=[0 1 0];
cMax=[1 0 1];

STEP=5;

maxUts(2)=0.1970;

for S=2 %[2 7 8];
    figure(S)
    clf
    hold on
    for k=1:length(c)
        r=regexp(c(k).name,['outputsS',num2str(S,'%1.2i')]);
        if ~isempty(r)
            st=open(['./org_data/',c(k).name]);

            X=st.Robot_arm.smooth_pos;
            V=st.Robot_arm.vel;
            t=st.Robot_arm.time;

            A=[gradient(V(:,1),t) gradient(V(:,2),t) gradient(V(:,3),t)];

            [Ut,Ux,Xgrid]=makeVPF(t,X,V,A);

            Ut=Ut-min(Ut); %Just a reference anyways
            mUt=max(Ut);
            nUt=Ut/maxUts(S);

            lnUt=log(nUt+.01);
            
            if mUt>maxUts(S)
                maxUts(S)=mUt;
            end
            
            for kk=STEP+1:length(t)
                plot3(X(kk-STEP:kk,1),X(kk-STEP:kk,2),X(kk-STEP:kk,3),'color',cMin*(1-lnUt(kk))+cMax*lnUt(kk))
            end
        end
    end
    
    axis equal
end
