function dissect(S)

global fJ

figure(575)
clf

qscale=1000;

S

NSUB=2;

subplot(NSUB,1,1)
hold on
plot(S.X(:,1),S.X(:,2),'b',S.X(:,1),S.Y(:,1),'g')
axis equal

Fe=zeros(size(S.X,1),2);
Fev=Fe;
Fff=Fe;
Fa=Fe;
Fi=Fe;
FB=Fe;
for k=1:size(S.X,1)
    qr=ikin(S.X(k,1:2));
    fJqrI=fJ(qr)';
    Fe(k,:)=(-fJqrI\S.E(k,1:2)')';
    Fev(k,:)=(-fJqrI\S.E(k,3:4)')';
    Fa(k,:)=(-fJqrI\S.oT(k,5:6)')';
    Fi(k,:)=(-fJqrI\S.oT(k,3:4)')';
    FB(k,:)=(-fJqrI\(S.Eb(k,1:2)'+S.Eb(k,3:4)'))';
    qd=ikin([S.X(k,1) S.Y(k,1)]);
    fJqdI=fJ(qd)';
    Fff(k,:)=-(-fJqdI\S.oT(k,1:2)')'; %Note, this is a sign quirk from solving for this most of the time
end
Fsum=Fff+Fi+Fa;

quiver(S.X(:,1),S.X(:,2),10*Fe(:,1)/qscale,10*Fe(:,2)/qscale,0,'Color','r')
quiver(S.X(:,1),S.X(:,2),Fev(:,1)/qscale,Fev(:,2)/qscale,0,'Color',[1 .5 .5])
quiver(S.X(:,1),S.X(:,2),Fa(:,1)/qscale,Fa(:,2)/qscale,0,'Color',[.5 .5 .5])
quiver(S.X(:,1),S.X(:,2),Fff(:,1)/qscale,Fff(:,2)/qscale,0,'Color','b')
quiver(S.X(:,1),S.X(:,2),Fi(:,1)/qscale,Fi(:,2)/qscale,0,'Color','g')
quiver(S.X(:,1),S.X(:,2),S.X(:,5)/qscale,S.X(:,6)/qscale,0,'Color',[0 .5 0])
quiver(S.X(:,1),S.X(:,2),Fsum(:,1)/qscale,Fsum(:,2)/qscale,0,'Color','m')
quiver(S.X(:,1),S.X(:,2),FB(:,1)/qscale,FB(:,2)/qscale,0,'Color','c')

plot(-.1*[1 1],.51+[0 10/qscale],'Color','k')

legend('Hand','Ergodic Desired','Error*10','Vel Error','Applied','Feedforward (Ergodic)','Inertial + Coriolis','Acceleration','Sum(ff + inertial + applied','Burdet FB')

subplot(NSUB,1,2)
vmFs=vecmag(Fsum);
vmFe=vecmag(Fe);
W=vmFe\vmFs
K=[Fe Fev]\Fsum
x=linspace(min(vmFe),max(vmFe),10);
plot(vmFe,vmFs,'b-.',x,W*x,'r')
ylabel('|Force Remainder|')
xlabel('|Error|')