function [summed,kerns]=supMJP(x0,dx,ti,tf,t)

summed=ones(length(t),1)*x0;
kerns=zeros(length(t),length(ti));

ldx=size(dx,2);

for k=1:size(dx,1)
    [kern,fop]=slmj5op(t,(tf(k)+ti(k))/2,tf(k)-ti(k));
    kerns(:,k)=kern;
    for der=1:3
        summed(:,ldx*(der-1)+1:ldx*der)=summed(:,ldx*(der-1)+1:ldx*der)+fop(:,der)*dx(k,:);
    end
end