function y=assembleLumps(t,lumps)

y=zeros(length(t),2);

for k=1:(length(lumps)/4)
    offset=4*(k-1);
    L=lumps(offset+(1:2))';
    C=lumps(offset+3);
    S=lumps(offset+4);
    
    tau=(t-C)/S+.5;
    tau=max(min(tau,1),0);
    kappa=(30*tau.^2-60*tau.^3+30*tau.^4);
    y=y+kappa*(L/S);
end
