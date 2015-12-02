%% Housekeeping

clc
clear all
close all

nbins=30;
datafile='free_exp_05stroke.mat';
%datafile='free_exp_MF.mat'

load(datafile)
en=dot(v',v');

%Plot the upside down bowl as a sanity check
figure(1) 
clf
inds=1:20:size(x,1);
plot3(x(inds,1),x(inds,2),en(inds),'.','markersize',.01)


%% 2D Histogram 

%Set up a grid
lower=min(x);
upper=max(x);
range=upper-lower;

selectors=[(x(:,1)-lower(1))/range(1), (x(:,2)-lower(2))/range(2)];
selectors=floor(selectors*nbins)+1;
for k=1:2
    selectors(selectors(:,k)==nbins+1,k)=nbins;
end

%Get the histogram
rawcounts=zeros(nbins);
fixedcounts=zeros(nbins);
for k=1:nbins
    f=find(selectors(:,1)==k);
    for kk=1:nbins
        f2=find(selectors(f,2)==kk);
        rawcounts(k,kk)=length(f2);
        velocities{k,kk}=v(f(f2),:);
    end
    counts(k,:)=rawcounts(k,:);
    counts(k,counts(k,:)==0)=1;
end

% This step also gets a sanity check
figure(2)
clf
surf(counts)

%% Deal with potentials

U=log(counts);
U=U-min(U(:));
U=U/max(U(:));

figure(3)
clf
surf(U)

%% Flatten velocity after a sanity check
totU=sum(U(:));
Uf=U(:);
ln=linspace(0,1,nbins);
[selmatX,selmatY]=meshgrid(ln*range(1)+lower(1),ln*range(2)+lower(2));
selmatXf=selmatX(:);
selmatYf=selmatY(:);
centroidX=sum(Uf.*selmatXf)/totU;
centroidY=sum(Uf.*selmatYf)/totU;
centroid=[centroidX,centroidY]
centroid=[.13 -.53]

velf=velocities(:);
distf=zeros(length(velf),1);
vel=distf;
for k=1:length(Uf)
    vel(k)=mean(vecmag(velf{k}).^2);
    distf(k)=norm([selmatXf(k) selmatYf(k)]-centroid)^2;
end
X=[vel distf ones(size(vel))];
[b,bint]=regress(Uf,X)
f=find(Uf>.7);
[b,bint]=regress(Uf(f),X(f,:))

figure(4)
clf
hold on
%plot(Uf,X*b,'r.')
plot(vel,Uf,'b.','markersize',.01)
%plot([0 1],[0 1],'m-')
axis equal

%% Different approach...show L^2 directly?
figure(5)
clf
f=find(Uf>3);
[b,bint]=regress(Uf(f),X(f,[1 3]))

velmat=zeros(nbins);
for k=1:nbins
    for kk=1:nbins
        velmat(k,kk)=mean(vecmag(velocities{k,kk}).^2);
        if isnan(velmat(k,kk))
            velmat(k,kk)=0;
        end
    end
end
surf(selmatX,selmatY,U+b(1)*velmat)
return