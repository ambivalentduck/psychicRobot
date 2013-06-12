clc
clear all

len=730;
signbalance=.5;
onevtwo=1/3;
y=zeros(len,1);
y(1)=0;

harshness=5; %Trade-off is structure vs balance. In theory, dynamic harshness dependent on k is the best trade-off and shortest methods section

%% vals "key":
%x
%y
%early pulse mag and +/- direction
%late pulse mag and +/- direction
%pseudorandom white noise magnitude
%shape (-1=none,0=triangle,1=square, 2=circle, 3=inf desired)
%cursor shown (0/1 false/true)

cost=inf;

while cost>3
    for k=2:len
        if y(k-1)==0
            if rand<signbalance
                y(k)=-1;
                signbalance=signbalance-harshness*1/len;
            else
                y(k)=1;
                signbalance=signbalance+harshness*1/len;
            end
        else
            if rand<onevtwo
                y(k)=0;
                onevtwo=onevtwo-harshness*2/len;
            else
                y(k)=-y(k-1);
                onevtwo=onevtwo+harshness*1/len;
            end
        end
    end
    diffs=y(2:end)-y(1:end-1);
    signbalance_=sum(diffs<0)-sum(diffs>0)
    onevtwo_=sum(abs(diffs)==1)-sum(abs(diffs)==2)
    cost=abs(onevtwo_)+abs(signbalance_);
end



x=[y [0;sign(diffs)] [0;abs(diffs)]];

N=5;
DIR=[-1 1];
DIST=[1 2];
ntypes=5;

%Make a deck
c=0;
for k=1:N
    for kk=1:length(DIR)
        for kkk=1:length(DIST)
            for kkkk=1:ntypes
                c=c+1;
                catme(c).dat=[DIR(kk); DIST(kkk); kkkk];
            end
        end
    end
end
deck=[catme.dat];
deck=deck(:,randperm(size(deck,2)));

c=0;
d=1;
k=inf;
sd2=size(deck,2);
dist=zeros(sd2,1);
while c<len
    c=c+1;
    k=k+1;
    if (k<5)||(c<30)||(d>sd2)
        out(c).dat=[y(c); .5; 0; 0; 0; -1; 1];
        continue
    end

    %if it's k>=7, reshuffle the deck and try again
    if (sum(x(c,2:3)'==deck(1:2,d))~=2)&&(k>5)
        deck(:,d:end)=deck(:,d-1+randperm(sd2-d+1));
    end
    if sum(x(c,2:3)'==deck(1:2,d))==2
        out(c).dat=[y(c); .5; 0; 0; 0; -1; 1];
        continue
    end

    %[c,k]
    dist(d)=c;
    k=0;
    switch deck(3,d)
        case 1
            out(c).dat=[y(c); .5; 1; 0; 0; -1; 1];
        case 2
            out(c).dat=[y(c); .5; -1; 0; 0; -1; 1];
        case 3
            out(c).dat=[y(c); .5; 0; 1; 0; -1; 1];
        case 4
            out(c).dat=[y(c); .5; 0; -1; 0; -1; 1];
        case 5
            out(c).dat=[y(c); .5; 0; 0; 1; -1; 1];
    end
    d=d+1;
end

figure(1)
clf
plot(1:len,y,dist,y(dist),'rx')
length(out)
size([out.dat])

[h]=hist(dist(2:end)-dist(1:end-1),5:11);
[(5:11)' h']

out(c).dat=[0; 0; 0; 0; 0; 0; 1]; %Extra copy of the first shape to allow target acquisition
for s=0:3 %5 warmups on each, cursor shown
    for k=1:5
        c=c+1;
        out(c).dat=[0; 0; 0; 0; 0; s; 1];
    end
end

for s=0:3 %5 warmups on each, cursor not shown
    for k=1:5
        c=c+1;
        out(c).dat=[0; 0; 0; 0; 0; s; 0];
    end
end

SnMag=zeros(2,4*5*20);
k=0;
for s=0:3
    for n=linspace(0,3,2.5)
        for ITER=1:20
            k=k+1;
            SnMag(:,k)=[n;s];
        end
    end
end
SnMag=SnMag(:,randperm(4*5*20));

for k=1:4*5*20
    c=c+1;
    out(c).dat=[0; 0; 0; 0; SnMag(:,k); 0];
end

o=[out.dat]';
o(1:730,1)=o(1:730,1)*.15-.025;
o(:,3:4)=o(:,3:4)*30;
o(1:730,5)=o(1:730,5)*2;
o=[(1:length(out))' o];
fid=fopen('../Data/input.dat','w');
fprintf(fid,'%5.0f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%1.0f\t%1.0f\n',o');
fclose(fid);
