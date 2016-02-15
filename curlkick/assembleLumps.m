function y=assembleLumps(t,lumps)

y=zeros(length(t),2);
for k=1:length(lumps)
    tau=(t-lumps(k).C)/lumps(k).S+.5;
    tau=max(min(tau,1),0);
    kappa=(30*tau.^2-60*tau.^3+30*tau.^4);
    y=y+kappa*(lumps(k).L/lumps(k).S);
end
