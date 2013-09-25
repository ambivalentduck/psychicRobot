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
%hash=[0 0 0 0 0]';
hash=char('A'+hash)';
name=['dynamicSimFile',hash]

fh=fopen([name,'.m'],'w');
fprintf(fh,['function simmed=',name,'(t,xvaf)\n\nglobal l1 l2 m1 m2 lc1 lc2 I1 I2 x0 kp0 kp1 kpgain real\n\ntic\n\nfigure(2)\nclf\nhold on\nplot(xvaf(:,1),xvaf(:,2),''b.'')\n\n']);

vary=[1.05 1.1 1.15 1.2 1.3];

count=0;
total=length(vary)*length(f)*2;
for v=1:length(vary)
    for k=1:length(f)
        for pm=1:2
            %Set gross parameters every time
            fprintf(fh,['l1=',num2str(l1),';\n']);
            fprintf(fh,['l2=',num2str(l2),';\n']);
            fprintf(fh,['mass=',num2str(mass),';\n']);
            fprintf(fh,['x0=[',num2str(x0(1)),' ',num2str(x0(2)),'];\n']);
            fprintf(fh,'kpgain=1;\nkp0=[10.8 2.83; 2.51 8.67];\nkp1=[3.18 2.15; 2.34 6.18];\n');
            if strcmp(simme.(f{k}).apply,'gross')&&(v==1) %Easiest just to overwrite whatever value was there and not think about which
                if pm==1
                    fprintf(fh,[simme.(f{k}).name,'=',num2str(simme.(f{k}).upper),';\n']);
                else
                    fprintf(fh,[simme.(f{k}).name,'=',num2str(simme.(f{k}).lower),';\n']);
                end
            end
            fprintf(fh,'lc1=.436*l1;\nlc2=.682*l2;\nm1=.028*mass;\nm2=.022*mass;\n');
            if strcmp(simme.(f{k}).apply,'fine')
                if pm==1
                    fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'*',num2str(vary(v)),';\n']);
                else
                    fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'/',num2str(vary(v)),';\n']);
                end
            end
            fprintf(fh,'I1=m1*(.322*l1)^2;\nI2=m2*(.468*l2)^2;\n');
            if strcmp(simme.(f{k}).apply,'end')
                if pm==1
                    fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'*',num2str(vary(v)),';\n']);
                else
                    fprintf(fh,[simme.(f{k}).name,'=',simme.(f{k}).name,'/',num2str(vary(v)),';\n']);
                end
            end
            fprintf(fh,['y=extract(t,xvaf,''reflex'');\nC=rand(1,3);\nplot(y(:,1),y(:,2),''Color'',C)\n']);
            fprintf(fh,['simmed(',num2str(k),',',num2str(v),',',num2str(pm),').y=y;\n']);
            count=count+1;
            fprintf(fh,['progress=',num2str(count/total),'\naxis equal\ndrawnow\ntoc\n\n']);
        end
    end
end
fclose(fh);

simmed=feval(name,t,xvaf);

for k=1:length(f)
    for v=1:length(vary)
        for pm=1:2
            simmed(k,v,pm).rms=sqrt(mean(vecmag2(yex-simmed(k,v,pm).y)))*1000;
        end
    end
end

fh=fopen('report.txt','w');

columns(1).v='Variable';
for k=1:length(vary)
    n=num2str(vary(k));
    columns(k+1).v=['\t/',n,'\t*',n];
end

fprintf(fh,[columns.v,'\n']);

rms=[[simmed(:,end,1).rms]' [simmed(:,end,2).rms]'];
v=max(rms,[],2);
[s,i]=sort(v);

for k=1:length(i)
    clear columns
    columns(1).v=f{i(k)};
    for v=1:length(vary)
        columns(v+1).v=['\t',num2str(simmed(i(k),v,1).rms),'\t',num2str(simmed(i(k),v,2).rms)];
    end
    fprintf(fh,[columns.v,'\n']);
end
fclose(fh);

%% Step 4: Monte Carlo Variance estimation

% First step, set up nominal values.

l1nom=.33;
l2nom=.34;
weightnom=175; %lbs
massnom=weightnom*.4535; %kg
shouldernom=[0 .45];

%Winters (1990)
lc1nom=.436*l1nom;
lc2nom=.682*l2nom;

m1nom=.028*massnom;
m2nom=.022*massnom;

%rog of gyration numbers from winters, rog=sqrt(I/m)
I1nom=m1nom*(.322*l1nom)^2;
I2nom=m2nom*(.468*l2nom)^2;

%Shoulder location
x0=origin+shouldernom; %Shoulder is measured in room coordinates relative to the workspace center

kpgain=1;
kp0nom=[10.8 2.83; 2.51 8.67];
kp1nom=[3.18 2.15; 2.34 6.18];

% Next, which parameters are really independent? Certainly not the
% terms inside a stiffness matrix...
% l1 l2 lc1 lc2 m1 m2 i1 i2 kp0 kp1... but really you want to fuzz up
% winters' numbers.

p=sobolset(20,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

varied=p(1:1000,:); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.
%The number after p controls the number of points. Remember, this N*10 is the number of individual sims done.

consistent=varied(:,1:10);
individual=varied(:,11:20);

range=.1*[l1 l2 lc1 lc2 m1 m2 I1 I2]; %Have to vary about nominal values to avoid correlated variance

sc1=size(consistent,1);

tic
for N=1:size(consistent,1)
    v=consistent(N,:);
    v(1:8)=(v(1:8)-.5).*range;
    v(9:10)=(v(9:10)-.5)*.4+1;
    l1=l1nom+v(1);
    l2=l2nom+v(2);
    lc1=lc1nom+v(3);
    lc2=lc2nom+v(4);
    m1=m1nom+v(5);
    m2=m2nom+v(6);
    I1=I1nom+v(7);
    I2=I2nom+v(8);
    kp0=kp0nom*v(9);
    kp1=kp1nom*v(10);
    simmedA(N).y=extract(t,xvaf,'reflex');

    v=individual(N,:);
    v(1:8)=(v(1:8)-.5).*range;
    v(9:10)=(v(9:10)-.5)*.4+1;
    l1=l1nom+v(1);
    l2=l2nom+v(2);
    lc1=lc1nom+v(3);
    lc2=lc2nom+v(4);
    m1=m1nom+v(5);
    m2=m2nom+v(6);
    I1=I1nom+v(7);
    I2=I2nom+v(8);
    kp0=kp0nom*v(9);
    kp1=kp1nom*v(10);
    simmedB(N).y=extract(t,xvaf,'reflex');

    for col=1:10
        v=consistent(N,:);
        v(col)=individual(N,col); %Overwrite just a single value at a time.
        v(1:8)=(v(1:8)-.5).*range;
        v(9:10)=(v(9:10)-.5)*.4+1;
        l1=l1nom+v(1);
        l2=l2nom+v(2);
        lc1=lc1nom+v(3);
        lc2=lc2nom+v(4);
        m1=m1nom+v(5);
        m2=m2nom+v(6);
        I1=I1nom+v(7);
        I2=I2nom+v(8);
        kp0=kp0nom*v(9);
        kp1=kp1nom*v(10);
        simmedAB(N,col).y=extract(t,xvaf,'reflex');
    end
    [N/sc1 toc/N ((sc1/N-1)*(toc))/60]
end
