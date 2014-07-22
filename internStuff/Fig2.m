clear;
clc;
load ./Data/output27.mat;
figure(1)
clf
hold on
%title('Subject 27')

for P=1:5
    for T=96*(P-1)+1:96*P
        if T<10
            continue
        end
        
x=trials(T).x;
y=trials(T).y;
c=trials(T).cursor;

xoff=.5*(trials(T).targcat==0);
yoff=1-.3*P;

plot(x(:,1)+xoff,x(:,2)+yoff,'k')
plot(y(:,1)+xoff,y(:,2)+yoff,'r')
plot(c(:,1)+xoff,c(:,2)+yoff,'linewidth',1.5)
    end
end

axis equal
axis off
 
text(-.3,1.2,'Baseline','HorizontalAlignment','center','VerticalAlignment','bottom')
text(-.3,.9,'Forces','HorizontalAlignment','center','VerticalAlignment','bottom')
text(-.3,.65,' Forces','HorizontalAlignment','center','VerticalAlignment','bottom')
text(-.3,.6,'+','HorizontalAlignment','center','VerticalAlignment','bottom')
text(-.3,.55,'Intent Shown','HorizontalAlignment','center','VerticalAlignment','bottom')
text(-.3,.3,'Forces','HorizontalAlignment','center','VerticalAlignment','bottom')
text(-.3,0,'Washout','HorizontalAlignment','center','VerticalAlignment','bottom')
text(0,1.4,'Outward','HorizontalAlignment','center','VerticalAlignment','bottom')
text(.5,1.4,'Inward','HorizontalAlignment','center','VerticalAlignment','bottom')


legend('Hand', 'Intent','Cursor')
legend('boxoff')

set(gcf,'color','w')