function out=genNewExp()

%% Step 0: Space parameters, experiment parameters
N=8;
D=3;
wdist=.1;
x=[-.21 .26];
y=[.37 .6];
early=[-1 1];
late=[-1 1];
white=[2];
angles=[60 180 300]*pi/180;

%% Step 1: Make a "deck" of perturbed reaches and shuffle it.
le=length(early);
ll=length(late);
lw=length(white);

lP=N*D*(le+ll+lw)
P=zeros(lP,4);
I=0;

for n=1:N
    for d=1:D
        for k=1:le
            I=I+1;
            P(I,:)=[early(k) 0 0 d];
        end
        for k=1:ll
            I=I+1;
            P(I,:)=[0 late(k) 0 d];
        end
        for k=1:lw
            I=I+1;
            P(I,:)=[0 0 white(k) d];
        end
    end
end
r=randperm(lP);
P=P(r,:);


%% Step 2: Take cards off the deck
p=1;
pr=0;
xy=[mean(x); mean(y)];

%vals "key":
%x
%y
%early pulse mag and +/- direction
%late pulse mag and +/- direction
%pseudorandom white noise magnitude
%shape (-1=none,0=triangle,1=square, 2=circle, 3=inf desired)
%direction metadata for easy deciphering later, -1 not "really" specified, 0=random, 1=60deg,%2=180deg, 3=300deg

dMemHack(1).v=[xy;0;0;0;-1;0];
%%grows like a simple array, not a matrix, trivial concatenation

while p<=30 %"warm up"
    p=p+1;
    [xy,d,dir]=walk(xy,angles,wdist,x,y);
    dMemHack(p).v=[xy;0;0;0;-1;dir];
end
pr=p;

lastexp=0;
deckind=1;
while deckind<=lP %Take cards off the stack as appropriate.
    p=p+1;
    lastexp=lastexp+1;
    [xy,d,dir]=walk(xy,angles,wdist,x,y);
    if (lastexp>5)&&(d>.01)&&(dir==P(deckind,4))
        dMemHack(p).v=[xy;P(deckind,1:3)';-1;dir];
        deckind=deckind+1;
        lastexp=0;
    else
        dMemHack(p).v=[xy;0;0;0;-1;dir];
    end
end
pr=p;


%% Step 3: Make a deck of shapes, shuffle it, and deal it.
for r=1:2
    for s=0:3
        for k=1:5
            pr=pr+1;
            dMemHack(pr).v=[0 0 0 0 0 s -1 ]';
        end
    end
end
SnMag=zeros(2,50*4);
for d=1:4
    for n=1:50
        SnMag(:,(d-1)*50+n)=[mod(n,10)/3;d];
    end
end
SnMag=SnMag(:,randperm(4*50));

for k=1:4*50
    dMemHack(pr+k).v=[0;0;0;0;SnMag(:,k);-1];
end

out=[dMemHack.v]';
out=[(1:size(out,1))' out]
fid=fopen('../Data/input.dat','w');
fprintf(fid,'%5.0f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%1.0f\t%1.0f\n',out');
fclose(fid);
    


function [p2,d,dir]=walk(p1,angles,wdist,x,y)
r=randperm(length(angles));
success=0;
for k=1:length(r)
    p2=p1+wdist*[cos(angles(r(k))); sin(angles(r(k)))];
    check=sum([p2(1)>x(1) p2(1)<x(2) p2(2)>y(1) p2(2)<y(2)]);
    if check==4
        d=minDist(x,y,p2);
        success=1;
        dir=r(k);
        break
    end
end
if ~success %if something goes very very wrong, recover gracefully
    p2=[x(1)+(x(2)-x(1))*rand;y(1)+(y(2)-y(1))*rand];
    d=minDist(x,y,p2);
    dir=0;
end

function dist=minDist(x,y,p)
d=zeros(4,1);
d(1)=pointlinedist(x,[y(1); y(1)],p);
d(2)=pointlinedist(x,[y(2); y(2)],p);
d(3)=pointlinedist([x(1); x(1)],y,p);
d(4)=pointlinedist([x(2); x(2)],y,p);
dist=min(d);

function dist=pointlinedist(x,y,p)
a=[x(1); y(1)];
n=[x(2); y(2)];
n=(n-a);
n=n/norm(n);
amp=a-p;
dist=norm(amp-dot(amp,n)*n);

