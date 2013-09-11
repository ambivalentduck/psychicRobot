function [l1, l2, shoulder,mass]=getSubjectParams(num)

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
    case {'3','9','10'}
        l1=.36;
        l2=.31;
        mass=225;
        shoulder=[0 .55];
    case {'11','12'}
        l1=.3;
        l2=.32;
        mass=213;
        shoulder=[0 .52];
end