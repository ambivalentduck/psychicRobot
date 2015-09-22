function [x1,x2,t1,t2]=makeSpaghetti

Ncand = 500;

xn = 1e-3;
subnum = 0;
t1(1) = 0;

t=0:.01:1;
xref=10*t.^3-15*t.^4+6*t.^5;
Uref=450*t.^4.*(t-1).^4;

fudge=10;

while (xn<.95)&&(subnum<20)
    subnum = subnum + 1;
    xnp1 = linspace(xn,1.05,Ncand)';
    
    Un = 506*exp(xn*-10.18)+1.1;
    Unp1 = 506*exp(xnp1.*-10.18)+1.1;
    
    %Vn=-450*xn^4*(xn - 1)^4;
    Vn=interp1(xref,Uref,xn,'nearest','extrap');
    %Vnp1=-450*xnp1.^4.*(xnp1 - 1).^4;
    Vnp1=interp1(xref,Uref,xnp1,'nearest','extrap');
    
    Anp1 = 589*exp(xnp1.*-9.37)+1.65;
    Tm2np1 = Unp1+expinv(rand(Ncand,1),Anp1);
    Tm2np1(Tm2np1> 0.05^-2)=0.05^-2;
    
    Rnp1 =fudge*(Vn-Vnp1)+((0.15.*(xnp1-xn)).^2.*Tm2np1)*1.875^2;

    pRnp1=exppdf(Rnp1,.2);
    pRnp1=pRnp1/sum(pRnp1);
    cpRnp1=cumsum(pRnp1);

    f=find(cpRnp1>=rand,1,'first');
    
    Tm2np1 = Tm2np1(f);
    t(subnum) = 1/sqrt(Tm2np1);
    x1(subnum) = 0.15*xn;
    xn = xnp1(f);
    x2(subnum) = 0.15*xn;
    t1(subnum+1) = t1(subnum) + t(subnum)/2;
    t2(subnum) = t1(subnum) + t(subnum);
end
t1(end)=[];

