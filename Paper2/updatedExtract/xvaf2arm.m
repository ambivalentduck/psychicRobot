function arm=xvaf2arm(xvaf)

global fJ getAlpha

arm=xvaf;

for k=1:size(xvaf,1)
    q=ikin(xvaf(k,1:2));
    fJq=fJ(q);
    qdot=fJq\xvaf(k,3:4)';
    qddot=getAlpha(q,qdot,xvaf(k,5:6)');
    torque=-fJq'*xvaf(k,7:8)';
    arm(k,:)=[q' qdot' qddot' torque'];
end

