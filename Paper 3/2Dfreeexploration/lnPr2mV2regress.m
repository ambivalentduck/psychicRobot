%% Housekeeping

clc
clear all
close all

nbins=20;
datafile='free_exp_05stroke.mat';
%datafile='free_exp_MF.mat'

load(datafile)
en=dot(v',v');
pow=dot(a',v');

%Plot the upside down bowl as a sanity check
figure(1) 
clf
inds=1:20:size(x,1);
plot3(x(inds,1),x(inds,2),pow(inds),'.','markersize',.01)


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
        pows{k,kk}=pow(f(f2));
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

figure(3)
clf
surf(U)

%% Flatten various values
Uf=U(:);
totU=sum(U(:));
ln=linspace(0,1,nbins);
[selmatX,selmatY]=meshgrid(ln*range(1)+lower(1),ln*range(2)+lower(2));
selmatXf=selmatX(:);
selmatYf=selmatY(:);
pows=pows(:);

meanpow=zeros(length(Uf),1);
for k=1:length(Uf)
    meanpow(k)=mean(abs(pows{k}));
end

figure(4)
clf
hold on
%plot(Uf,X*b,'r.')
plot(Uf,meanpow,'b.','markersize',.01)
%plot([0 1],[0 1],'m-')
axis equal
return
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