clc
clear all

% establish a straight line reach potential field

% R needs U
% R needs T
% Relationship between L and U
Ncand = 500;
figure(438)
clf
hold on

xn = 1e-3;
subnum = 0;
% t1 = zeros(subnum+1,1);
% t2 = zeros(subnum,1);
% x1 = zeros(subnum,1);
% x2 = zeros(subnum,1);
% t = zeros(subnum,1);
t1(1) = 0;
fudge=-1;
while (xn<1)&&(subnum<20)
    subnum = subnum + 1;
    xnp1 = linspace(xn,1.03,Ncand)';
    
    Un = 506*exp(xn*-10.18)+1.1;
    Unp1 = 506*exp(xnp1.*-10.18)+1.1;
    
    Vn=-450*xn^4*(xn - 1)^4;
    Vnp1=-450*xnp1.^4.*(xnp1 - 1).^4;
    
    Anp1 = 589*exp(xnp1.*-9.37)+1.65;
    Tm2np1 = Unp1+expinv(rand(Ncand,1),Anp1);
    Tm2np1(Tm2np1> 0.05^-2)=0.05^-2;
    
    Rnp1 = fudge*(Vnp1-Vn)+((0.15.*(xnp1-xn)).^2.*Tm2np1)*1.875^2;
    % figure(1)
    % clf
    % subplot(2,1,1);
    % plot(xnp1,1./sqrt(Tm2np1),'k')
    % subplot(2,1,2);
    % plot(xnp1,Rnp1,'b')
    % hold on
    % choose R from an expoential distribution temp = 0.67
    pRnp1=exppdf(Rnp1,.2);
    pRnp1=pRnp1/sum(pRnp1);
    cpRnp1=cumsum(pRnp1);
    % for i = 1:Ncand
    %     f(i)=find(cpRnp1>=rand,1,'first');
    % end
    f=find(cpRnp1>=rand,1,'first');
    % plot(xnp1(f),zeros(size(f))+0.2+rand(size(f))*0.2,'k.','MarkerSize',1)
    % plot(xnp1,cpRnp1,'m-')
    % plot(xnp1,gradient(cpRnp1)./gradient(xnp1),'g-')
    % xn = xnp1(f);
    plot([xn,xnp1(f)],subnum+rand*0.5+[0 0],'k-')
    Tm2np1 = Tm2np1(f);
    t(subnum) = 1/sqrt(Tm2np1);
    x1(subnum) = 0.15*xn;
    %xn = mean([xnp1(f),xn]);
    xn = xnp1(f);
    x2(subnum) = 0.15*xn;
    t1(subnum+1) = t1(subnum) + t(subnum)/2;
    t2(subnum) = t1(subnum) + t(subnum);
end
t1(end)=[];
plotSubefforts(t1,t2,x1,x2);
