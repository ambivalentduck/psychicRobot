function y=cosLumpy(t,lump)

if size(t,1)<size(t,2)
    t=t';
end

tau=(t-lump.C)/lump.S+.5;
tau=max(min(tau,1),0);
kappa=pi*sin(pi*tau)/2;
y=kappa*(lump.L/lump.S);