function y=lumpy(t,lump)

if size(t,1)<size(t,2)
    t=t';
end

tau=(t-lump.C)/lump.S+.5;
tau=max(min(tau,1),0);
kappa=(30*tau.^2-60*tau.^3+30*tau.^4);
y=kappa*(lump.L/lump.S);