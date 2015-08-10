function [fitparams,nominal]=fitAllParams(S)

global fJ getAlpha

load(['../Data/Data_pulse/pulse',num2str(S),'W.mat'])

set2dGlobals(params.l1, params.l2, params.origin, params.shoulder, params.mass)

%% Discard the empties
ss={storeme.starts};
notempty=zeros(length(ss),1);
for k=1:length(ss)
    notempty(k)=~isempty(storeme(k).starts);
end

storeme=storeme(notempty==1);

%% Compute q, qd, and tau_ext

subinds=41:80;

for SM=1:length(storeme)
    X=storeme(SM).X(subinds,:);
    Y=storeme(SM).Y(subinds,:);
    storeme(SM).qreal=zeros(6,size(X,1));
    storeme(SM).qdes=zeros(6,size(X,1));
    storeme(SM).tau_app=zeros(2,size(X,1));

    for k=1:size(X,1)
        %Qreal and torque
        q=ikin(X(k,1:2));
        fJq=fJ(q);
        qdot=fJq\X(k,3:4)';
        qddot=getAlpha(q,qdot,X(k,5:6)');
        storeme(SM).tau_app(:,k)=-fJq'*X(k,7:8)'; %Negative sign is consistent with extract.m
        storeme(SM).qreal(:,k)=[q; qdot; qddot];

        %Qdes
        q=ikin(Y(k,1:2));
        fJq=fJ(q);
        qdot=fJq\Y(k,3:4)';
        qddot=getAlpha(q,qdot,Y(k,5:6)');
        storeme(SM).qdes(:,k)=[q; qdot; qddot];
    end
end

%% Since we're ignoring time shifts, just concatenate

qr=[storeme.qreal];
qd=[storeme.qdes];
tau_a=[storeme.tau_app];

%% Plug everything into an optimizer
l1=params.l1;
l2=params.l2;
mass=params.mass;
m1=.028*mass;
m2=.022*mass;

kp0=[10.8 2.83; 2.51 8.67];
kp1=[3.18 2.15; 2.34 6.18];

nominal=[.436*l1,.682*l2,m1,m2,m1*(.322*l1)^2,m2*(.468*l2)^2,kp0(:)',kp1(:)'];
fitparams=deoptImpedanceParams(qd,qr,tau_a,l1,nominal);

