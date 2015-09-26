function [x1,x2,t1,t2]=makeSpaghetti

Ncand = 500;

xn = 1e-3;
subnum = 0;
t1(1) = 0;

t=0:.01:1;
xref=10*t.^3-15*t.^4+6*t.^5;
Uref=1*450*t.^4.*(t-1).^4;

targAccuracy=.01;
reachL=.15;
massfudge=.5;

while (xn<(1-targAccuracy))&&(subnum<20)
    subnum = subnum + 1;
    xnp1 = linspace(xn,1+targAccuracy,Ncand)';
    
    Un = 506*exp(-10.18*xn)+1.1;
    Unp1 = 506*exp(-10.18*xnp1)+1.1;
    
    Vn=interp1(xref,Uref,xn,'nearest','extrap');
    Vnp1=interp1(xref,Uref,xnp1,'nearest','extrap');
    
    Anp1 = 589*exp(-9.37*xnp1)+1.65;
    Tm2np1 = Unp1+expinv(rand(Ncand,1),Anp1);
    %Tm2np1(Tm2np1> 0.05^-2)=0.05^-2;
    
    Rnp1 = (Vn-Vnp1)+massfudge*((reachL.*(xnp1-xn)).^2.*Tm2np1)*1.875^2;

    pRnp1=exppdf(Rnp1,.6);
    pRnp1=pRnp1/sum(pRnp1);
    cpRnp1=cumsum(pRnp1);

    f=find(cpRnp1>=rand,1,'first');
    
    t(subnum) = 1/sqrt(Tm2np1(f));
    x1(subnum) = reachL*xn;
    xn = xnp1(f);
    x2(subnum) = reachL*xn;
    t1(subnum+1) = t1(subnum) + t(subnum)/2;
    t2(subnum) = t1(subnum) + t(subnum);
end
t1(end)=[];

