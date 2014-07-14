clc
clear all

syms t t1 t2 t3 t4 real

ta=(t-t1)/(t2-t1)
k1=10*ta.^3-15*ta.^4+6*ta.^5

tb=(t-t3)/(t4-t3)
k2=10*tb.^3-15*tb.^4+6*tb.^5

%Unless t2>t3, these things don't touch so multiplication = 0

ip=simplify(int(k1*k2,t,t3,t2))

ipv=inline(vectorize(ip))

%Choose t1=0 because that's just a shift on t.
%Choose t4=1 because that's just division by t4
%Now it's 3 dimensions: t2, t3, ip

area=0:.01:1;

[X,Y]=meshgrid(area,area);

Z=ipv(0,X,Y,1);
Z(Y>=X)=0;

surf(X,Y,Z)
