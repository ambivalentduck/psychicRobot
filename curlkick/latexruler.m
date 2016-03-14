clc
clear all

width=8; %3.4*2;

f=figure(1);
clf
set(f,'units','inches');
set(f,'Position',[1,1,width,1])
plot(0:width,ones(width+1,1),'kx')
for k=0:width+1
    text(k,1,num2str(k))
end
set(gca,'position',[0 0 1 1])