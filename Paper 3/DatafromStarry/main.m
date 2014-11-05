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

for S=subs
    stringS=num2str(S,'%1.2i');
    for k=1:3
        pre(k)=load(['../org_data/outputsS',stringS,'D1F',num2str(k),'_org_session.mat']);
        post(k)=load(['outputsS',stringS,'D6F',num2str(k),'_org_session.mat']);
        
        
    end
    
end