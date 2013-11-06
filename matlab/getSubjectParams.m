function params=getSubjectParams(num)

switch num
    case {'1', '500', '21', '31','32','55','729'} %me
        l1=.33;
        l2=.34;
        mass=175;
        shoulder=[0 .45];
    case {'300'} %Yazan
        l1=.31;
        l2=.34;
        mass=210;
        shoulder=[0 .51];
    case {'301'} %Guy
        l1=.33;
        l2=.345;
        mass=183;
        shoulder=[0 .49];
    case {'3','9','10'} %Jim
        l1=.36;
        l2=.31;
        mass=225;
        shoulder=[0 .55];
    case {'11','12'} %Tes
        l1=.3;
        l2=.32;
        mass=213;
        shoulder=[0 .52];
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