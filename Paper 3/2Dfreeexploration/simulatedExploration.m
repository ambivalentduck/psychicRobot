clc
clear all

%Syntax is maxV=saccadefunc; [x,v]=saccadefunc(t0,tf,t);
%saccadefunc=@minjerkSaccade; %Closer to reality
saccadefunc=@linearSaccade; %Important for testing since E(t) is constant
mass=1;
temp=.01;
bounty=0;
nbins=10;
timetemp=.2;

sacMaxV=saccadefunc();
sacRu=.5*mass*sacMaxV^2;

%% Set up a 2D potential field

[X,Y]=meshgrid(linspace(-.2,.2,nbins),linspace(.35,.7,nbins)); %Mimic manipulandum

U=(X-0).^2+(Y-.5).^2; %The most important decision in this file
%U=0*U;

figure(1)
clf
surf(X,Y,U)

%% Explore the space with jumps

jumps(1).tf=0;
jumps(1).xf=[0 .5];
jumps(1).U=0; %U = sum deltaU, so avoid doing lookups.
k=1;


while k<10000 %jumps(k).tf<10 %Seconds to run the simulation
    %Each point on the grid is a candidate, deal in unpixelated reality later
    
    L2=(X-jumps(k).xf(1)).^2+(Y-jumps(k).xf(2)).^2;
    deltaU=U-jumps(k).U;
    
    % Tn2=(bounty-deltaU)./(3*L2*sacRu); %Based on Shadmehrian T_crit
    Tn2=0+exprnd(timetemp,size(L2));
    J=exp(-(deltaU-bounty+sacRu*L2.*Tn2)/temp);
    %J(J<0)=0;
    %J((Tn2<=0)|(L2==0))=0; %Needed for T_crit approach
    %Choose among the candidates with Pr prop to e^-J/temp
    cJ=cumsum(J(:));
    i=find(cJ>rand*cJ(end),1,'first');
    
    %Store the jump
    jumps(k+1).x0=jumps(k).xf;
    jumps(k+1).xf=[X(i) Y(i)];
    jumps(k+1).L2=L2(i);
    jumps(k+1).Tn2=Tn2(i);
    jumps(k+1).deltaU=deltaU(i);
    jumps(k+1).U=jumps(k).U+deltaU(i);
    
    k=k+1;
end

%% Do an immediate sanity check

figure(2)
clf

xfs=vertcat(jumps(:).xf)+.002*randn(k,2);

plot(xfs(:,1),xfs(:,2),'.')
axis equal

%% Ask the critical question: is Tn2 exponential?

figure(3)
clf
subplot(2,1,1)
tn2=[jumps.Tn2];
[f,x,lo,hi]=ecdf(tn2);
hold on
h=fill([x(~isnan(lo)); wrev(x(~isnan(hi)))],[lo(~isnan(lo)); wrev(hi(~isnan(hi)))],'k');
plot(x,f,'k')
plot(x,expcdf(x,expfit(tn2)),'r')

subplot(2,1,2)
en=[jumps.L2].*[jumps.Tn2];
[f,x,lo,hi]=ecdf(en);
hold on
h=fill([x(~isnan(lo)); wrev(x(~isnan(hi)))],[lo(~isnan(lo)); wrev(hi(~isnan(hi)))],'k');
plot(x,f,'k')
plot(x,expcdf(x,expfit(en)),'r')