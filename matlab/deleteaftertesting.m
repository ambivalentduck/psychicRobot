
clc
clear all

figure(1)
clf
hold on

X=(1:40)';

Y=logspace(-2, 2, 100)

for k=Y
    a=getR2(X,X+k*randn(40,1));
    plot(k,a(1),'r.')
    plot(k,a(2),'b.')
end