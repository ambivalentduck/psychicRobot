function params=getSubjectParams(num)

switch num
    case '300'
        l1=.31;
        l2=.34;
        mass=210;
        shoulder=[0 .51];
        origin=[-0.0095, 0.4250];
    case '301'
        l1=.33;
        l2=.345;
        mass=183;
        shoulder=[0 .49];
        origin=[-0.0095, 0.4250];
    case '324'
        l1=.28;
        l2=.32;
        shoulder=[0 .53];
        mass=213;
        origin=[-0.0095, 0.4250];
    case '789'
        l1=.28;
        l2=.3;
        shoulder=[-.03 .48];
        mass=120;
        origin=[-0.0095, 0.4250];
    case '5'
        l1=.34;
        l2=.34;
        shoulder=[-0.02, .55];
        mass=208;
        origin=[0, .44];
    case '6'
        l1=.37;
        l2=.38;
        shoulder=[-0.03, .51];
        mass=190;
        origin=[0, .44];
    case '7'
        l1=.28;
        l2=.31;
        shoulder=[-.05, .48];
        mass=140;
        origin=[0, .44];
    case '8'
        l1=.29;
        l2=.33;
        shoulder=[0, .5];
        mass=160;
        origin=[0, .44];
end

params.l1=l1;
params.l2=l2;
params.shoulder=shoulder;
params.origin=origin;

params.mass=0.453592*mass; %Kg vs lbs
