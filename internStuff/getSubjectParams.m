function params=getSubjectParams(num)

switch num
    case {'1'} %me
        l1=.33;
        l2=.34;
        mass=175;
        shoulder=[.025 .92-.4378];
    case {'2'} %tejas
        l1=.34;
        l2=.36;
        mass=145;
        shoulder=[0 .99-.4378];
    case {'3'} %christine
        l1=.33;
        l2=.34;
        mass=120;
        shoulder=[.14 .96-.4378];
end

params.l1=l1;
params.l2=l2;
params.shoulder=shoulder;
params.origin=[0, 0.4378];

if exist('mass','var')
    params.mass=0.453592*mass;
else
    params.mass=((l1+l2)/.77)*91; %Shitty locally linear scaling of weight with height, remember that 200 lbs = 91 kg.
end