clc
clear all

figure(26)
clf
hold on
plot(-1*[1 1],-1*[1 1],'b-')
plot(-1*[1 1],-1*[1 1],'r-')

E=0:.01:1;

for n=[1 1.5 2 3]
    PE=gampdf(E,n,.15);
    plot(E,PE)
    text(E(40),PE(40),['n=',num2str(n)],'horizontalalignment','left')
end

for n=[1 1.5 2 3]
    PE=gampdf(E,n,.5);
    plot(E,PE,'r')
    text(E(20),PE(20),['n=',num2str(n)],'horizontalalignment','right')
end

xlim([0 1])
ylim([0 3.5])
xlabel('Energy, J')
ylabel('Pr(Energy)')
legend('T=.35','T=.7')