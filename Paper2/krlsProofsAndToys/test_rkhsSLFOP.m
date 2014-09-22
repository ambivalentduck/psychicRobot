clc
clear all

%Step 0 to everything else is knowing that support-limited 5th order
%polynomials are a RKHS.

%Rule 1: symmetry; h*g=g*h

tc=0;
ts=1;

testvec=-1.5:.5:1;
testsuccess=zeros(length(testvec));

for ht=1:length(testvec)
    h=slfop(testvec(ht),tc,ts);
    for gt=1:length(testvec)
        g=slfop(testvec(gt),tc,ts);
        testsuccess(ht,gt)=h*g-g*h;
    end
end

these_should_all_be_zero=testsuccess

if sum(testsuccess(:))==0
    disp('Symmetry Test passed.')
end

%Rule 2: scaling and distributive property; (cf+dg)*h=c(f*h)+d(g*h)

testsuccess=zeros(length(testvec),length(testvec),length(testvec));

syms c d real

for ht=1:length(testvec)
    h=slfop(testvec(ht),tc,ts);
    for gt=1:length(testvec)
        g=slfop(testvec(gt),tc,ts);
        for ft=1:length(testvec)
            f=slfop(testvec(ft),tc,ts);
            testsuccess(ht,gt,ft)=(c*f+d*g)*h~=c*(f*h)+d*(g*h);
        end
    end
end

these_should_all_be_zero=max(testsuccess,[],3)

if sum(testsuccess(:))==0
    disp('Scaling and Distributive Test passed.')
end