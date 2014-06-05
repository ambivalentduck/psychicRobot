function params=getSubjectParamsCurl(num)

switch num
    case '1'
        l1=.245;
        l2=.335;
        shoulder=[0 .45];
        mass=210;
    case '2'
        l1=.255;
        l2=.33;
        shoulder=[0 .44];
        mass=150;
    case '3'
        l1=.25;
        l2=.31;
        shoulder=[0 .37];
        mass=250;
    case '4'
        l1=.265;
        l2=.335;
        shoulder=[0 .46];
        mass=170;
    case '5'
        l1=.295;
        l2=.335;
        shoulder=[0 .48];
        mass=145;
    case '6'
        l1=.275;
        l2=.315;
        shoulder=[.02 .45];
        mass=145;
    case '7'
        l1=.255;
        l2=.305;
        shoulder=[0 .43];
        mass=150;
    case '8'
        l1=.34;
        l2=.31;
        shoulder=[0 .45];
end

params.l1=l1;
params.l2=l2;
params.shoulder=shoulder;
params.origin=[-0.0095, 0.4250];

if exist('mass','var')
    params.mass=0.453592*mass;
else
    params.mass=((l1+l2)/.77)*91; %Shitty locally linear scaling of weight with height, remember that 200 lbs = 91 kg.
end

