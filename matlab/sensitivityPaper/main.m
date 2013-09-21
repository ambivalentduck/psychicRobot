clc
clear all

global measuredVals measuredTime l1 l2 m1 m2 lc1 lc2 I1 I2 x0 kp0 kp1 kpgain

%% Step 0: Set physical parameters

LEFT=.281;
RIGHT=-.3951;
TOP=.18;
BOTTOM=.68;
origin=[(LEFT+RIGHT)/2,(TOP+BOTTOM)/2]; %Might as well duplicate the physical workspace...why not?

%Use my own physical parameters because...why not? Also, no IRB/justification
l1=.33;
l2=.34;
weight=175; %lbs
mass=weight*.4535; %kg
shoulder=[0 .45];

%Winters (1990)
lc1=.436*l1;
lc2=.682*l2;

m1=.028*mass;
m2=.022*mass;

%rog of gyration numbers from winters, rog=sqrt(I/m)
I1=m1*(.322*l1)^2;
I2=m2*(.468*l2)^2;

%Shoulder location
x0=origin+shoulder; %Shoulder is measured in room coordinates relative to the workspace center

kp0=[10.8 2.83; 2.51 8.67];
kp1=[3.18 2.15; 2.34 6.18];

%% Step 1: Set up a basic kicked movement and forward simulate.

t=0:.005:2;
coeff=calcminjerk([0 .5],[.15 .5],[0 0],[0 0],[0 0],[0 0],0,.7);
[x,v,a]=minjerk(coeff,t);
x=x';
v=v';
a=a';

i=find(t>=.7);
x(i,1)=x(i(1)-1,1);
x(i,2)=x(i(1)-1,2);
v(i,:)=0;
a(i,:)=0;

f=zeros(length(t),2);

i=find((t>=.1)&(t<=.15));
f(i,2)=15;

xvaf=[x v a f];

measuredVals=xvaf;
measuredTime=t;

for k=1:size(xvaf,1)
    q=ikin(xvaf(k,1:2));
    fJq=fJ(q,l1,l2);
    qdot=fJq\xvaf(k,3:4)';
    qddot=getAlpha(q,qdot,xvaf(k,5:6)',l1,l2);
    torque=-fJq'*xvaf(k,7:8)';
    measuredVals(k,:)=[q' qdot' qddot' torque'];
end

q0=measuredVals(1,1:4);

%forward simulate
kpgain=1;
[T,X]=forwardReflexHelper(t,q0);
[T,XSM]=ode45(@armdynamicsShadMuss,t,q0);

y=zeros(length(T),6);

for k=1:length(T)
    y(k,1:2)=fkin(X(k,1:2));
    y(k,3:4)=(fJ(X(k,1:2),l1,l2)*X(k,3:4)')';
end

gT=gradient(t)';
y(:,5)=gradient(y(:,3))./gT;
y(:,6)=gradient(y(:,4))./gT;

figure(1)
clf
hold on
plot(x(:,1),x(:,2),'b',y(:,1),y(:,2),'k.')
quiver(y(:,1),y(:,2),f(:,1),f(:,2),'b')
axis equal

%% Step 2: Extract to demonstrate accuracy

yex=extract(t,[y f],'reflex');
plot(yex(:,1),yex(:,2),'r.')

legend('Intent','Forward Sim','Forces','Extracted Intent')

xvaf=[y f];

%% Step 3: Set up data for mass simulation

clear simme fineparams endparams grossparams

% Make the structure: gross params first
fieldname='l1';
simme.(fieldname).name='l1';
simme.(fieldname).plotName='l1';
simme.(fieldname).lower=l1-.03;
simme.(fieldname).upper=l1+.03;
simme.(fieldname).apply='gross';

fieldname='l2';
simme.(fieldname).name='l2';
simme.(fieldname).plotName='l2';
simme.(fieldname).lower=l2-.03;
simme.(fieldname).upper=l2+.03;
simme.(fieldname).apply='gross';

fieldname='x0_1';
simme.(fieldname).name='x0(1)';
simme.(fieldname).plotName='x0_1';
simme.(fieldname).lower=x0(1)-.03;
simme.(fieldname).upper=x0(1)+.03;
simme.(fieldname).apply='gross';

fieldname='x0_2';
simme.(fieldname).name='x0(2)';
simme.(fieldname).plotName='x0_2';
simme.(fieldname).lower=x0(2)-.03;
simme.(fieldname).upper=x0(2)+.03;
simme.(fieldname).apply='gross';

fieldname='mass';
simme.(fieldname).name='mass';
simme.(fieldname).plotName='mass';
simme.(fieldname).lower=mass-9; %20 pounds
simme.(fieldname).upper=mass+9;
simme.(fieldname).apply='gross';

% Make the structure, fine params in a batch of +/- 20% because Winters
% could be off...say...just in the forearm
fineparams={'lc1','lc2','m1','m2','kpgain'};
for k=1:length(fineparams)
    simme.(fineparams{k}).name=fineparams{k};
    simme.(fineparams{k}).plotName=fineparams{k};
    simme.(fineparams{k}).apply='fine';
end

endparams={'I1','I2'};
endparamnames={'I1','I2'};
for k=1:4
    endparams{2+k}=['kp0(',num2str(mod(k-1,2)+1),',',num2str(ceil(k/2)),')'];
    endparams{6+k}=['kp1(',num2str(mod(k-1,2)+1),',',num2str(ceil(k/2)),')'];
    endparamnames{2+k}=['kp0_',num2str(mod(k-1,2)+1),num2str(ceil(k/2))];
    endparamnames{6+k}=['kp1_',num2str(mod(k-1,2)+1),num2str(ceil(k/2))];
end
for k=1:length(endparams)
    simme.(endparamnames{k}).name=endparams{k};
    simme.(endparamnames{k}).plotName=endparamnames{k};
    simme.(endparamnames{k}).apply='end';
end


f=fieldnames(simme);

hash=floor(rand(5,1)*24+1);
hash=char('A'+hash)';
name=['dynamicSimFile',hash]

fh=fopen([name,'.m'],'w');
fprintf(fh,['function simmed=',name,'(t,xvaf)\n\nglobal l1 l2 m1 m2 lc1 lc2 I1 I2 x0 kp0 kp1 kpgain real\n\ntic\n\nfigure(2)\nclf\nhold on\nplot(xvaf(:,1),xvaf(:,2),''b.'')\n\n']);

VARY=.3;

for k=1:length(f)
    for v=1:2
        %Set gross parameters every time
        fprintf(fh,['l1=',num2str(l1),';\n']);
        fprintf(fh,['l2=',num2str(l2),';\n']);
        fprintf(fh,['mass=',num2str(mass),';\n']);
        fprintf(fh,['x0=[',num2str(x0(1)),' ',num2str(x0(2)),'];\n']);
        fprintf(fh,'kpgain=1;\nkp0=[10.8 2.83; 2.51 8.67];\nkp1=[3.18 2.15; 2.34 6.18];\n');
        if strcmp(simme.(f{k}).apply,'gross') %Easiest just to overwrite whatever value was there and not think about which
            if v==1
                fprintf(fh,[simme.(f{k}).name,'=',num2str(simme.(f{k}).upper),';\n']);
            else
                fprintf(fh,[simme.(f{k}).name,'=',num2str(simme.(f{k}).lower),';\n']);
            end
        end
        fprintf(fh,'lc1=.436*l1;\nlc2=.682*l2;\nm1=.028*mass;\nm2=.022*mass;\n');
        if strcmp(simme.(f{k}).apply,'fine')
            if v==1
                fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'*',num2str(1-VARY),';\n']);
            else
                fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'*',num2str(1+VARY),';\n']);
            end
        end
        fprintf(fh,'I1=m1*(.322*l1)^2;\nI2=m2*(.468*l2)^2;\n');
        if strcmp(simme.(f{k}).apply,'end')
            if v==1
                fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'*',num2str(1-VARY),';\n']);
            else
                fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'*',num2str(1-VARY),';\n']);
            end
        end
        fprintf(fh,['y=extract(t,xvaf,''reflex'');\nC=rand(1,3);\nplot(y(:,1),y(:,2),''Color'',C)\nsimmed(',num2str((2*(k-1)+v)),').y=y;\nsimmed(end).name=''',simme.(f{k}).plotName,''';\n']);
        fprintf(fh,['progress=',num2str((2*(k-1)+v)/(length(f)*2)),'\naxis equal\ndrawnow\ntoc\n\n']);
    end
end
fclose(fh);

simmed=feval(name,t,xvaf);



