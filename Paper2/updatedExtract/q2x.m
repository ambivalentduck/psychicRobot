function x=q2x(q)

global fJ

x=q;
for k=1:size(q,1)
    x(k,1:2)=fkin(q(k,1:2));
    fJq=fJ(q(k,1:2));
    x(k,3:4)=(fJq*q(k,3:4)')';
end