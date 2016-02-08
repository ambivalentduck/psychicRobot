function [l1,l2,x0]=getSubjectParams(k)

switch k
    case 1
        l1=.255;
        l2=.33;
        shoulder=[0 .46];
    case 2
        l1=.25;
        l2=.31;
        shoulder=[0 .37]; %???
    case 3
        l1=.265;
        l2=.335;
        shoulder=[0 .46];
    case 4
        l1=.295;
        l2=.335;
        shoulder=[0 .48];
    case 5
        l1=.275;
        l2=.315;
        shoulder=[.02 .45];
    case 6
        l1=.255;
        l2=.305;
        shoulder=[0 .43];
    case 7
        l1=.34;
        l2=.31;
        shoulder=[0 .45];
    case 8
        l1=.29;
        l2=.36;
        shoulder=[.027 .405];
end

LEFT=.325;
RIGHT=-.344;
TOP=.17;
BOTTOM=.68;
origin=[(LEFT+RIGHT)/2,(TOP+BOTTOM)/2];

x0=origin+shoulder;