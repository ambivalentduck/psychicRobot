function x=q2x(q)

for k=1:size(q,1)
    x(k,:)=fkin(q(k,1:2));
end