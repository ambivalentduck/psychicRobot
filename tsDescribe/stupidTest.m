function [h,p]=stupidTest(x)

x=x-mean(x);
x=x/std(x);

[h,p]=kstest(x);