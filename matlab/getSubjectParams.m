function [l1, l2, shoulder,mass]=getSubjectParams(num)

switch num
    case '1'
        l1=.33;
        l2=.34;
        mass=160;
        shoulder=[0 .54];
    case {'3','9','10'}
        l1=.36;
        l2=.31;
        mass=225;
        shoulder=[0 .55];
end