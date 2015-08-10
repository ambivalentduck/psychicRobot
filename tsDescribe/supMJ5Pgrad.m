function [cost, grad]=supMJ5Pgrad(P)

global xdot tfit

%% Step 0 P->w,tc,ts
datadim=min(size(xdot));

lPD=length(P)/(datadim+2);

r=reshape(P,datadim+2,lPD);

w=r(1:2,:)';
tc=r(3,:);
ts=abs(r(4,:));

%% Step 1 cost....sum error.^2 over tc(1)+/-ts(1)

uppers=tc+ts/2;
lowers=tc-ts/2;

inds=find((tfit>=min(lowers))&(tfit<=max(uppers)));
tinds=tfit(inds);

[ydot,kerns]=supMJ5P(w,tc,ts,tinds);

fit_error=ydot-xdot(inds,:);
gradw=zeros(lPD,datadim);
gradtc=zeros(lPD,1);
gradts=zeros(lPD,1);

for k=1:lPD
    ta=(tinds-tc(k))/ts(k)+.5;
    tcdiff=-(60*ta-180*ta.^2+120*ta.^3)/(ts(k)^2);
    
    td=(tinds-tc(k))/(ts(k)^2);
    tsdiff=-(kerns(:,k)+td.*(60*ta-180*ta.^2+120*ta.^3))/ts(k);
   
    for kk=1:datadim
        gradw(k,kk)=sum(fit_error(:,kk).*kerns(:,k));
        gradtc(k)=gradtc(k)+sum(w(k,kk)*fit_error(:,kk).*tcdiff);
        gradts(k)=gradts(k)+sum(w(k,kk)*fit_error(:,kk).*tsdiff);
    end
end

cost=sum(sum(fit_error.^2))/2;

grad=zeros(datadim+2,lPD);
for k=1:lPD
    grad(:,k)=[gradw(k,:),gradtc(k),gradts(k)]';
end
grad=grad(:);