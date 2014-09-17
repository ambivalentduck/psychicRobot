function [cost, grad]=sparseupMJPfit(P)

global xdot t

%% Step 0 P->w,ts

r=reshape(P,3,5);
w=r(:,1:3);
tc=r(:,4);
ts=r(:,5);


%% Step 1 cost....sum error.^2 over tc(1)+/-ts(1)

uppers=tc+ts/2;
lowers=tc-ts/2;

inds=find((t>lowers(1))&(t<uppers(1)));
tinds=t(inds);

ydot=zeros(length(inds),3);
kerns=zeros(length(inds),3);

for k=1:3
    ta=(tinds-tc(k))/ts+.5;
    kerns(:,k)=(30*ta.^2-60*ta.^3+30*ta.^4)/ts(k);
    kerns(tinds<=lowers(k),k)=0;
    kerns(tinds>=uppers(k),k)=0;
    ydot=ydot+kerns(:,k)*w(k,:);
end

overlap=max(0, min(uppers(2),uppers(3))-max(lowers(2),lowers(3)));

fit_error=xdot(inds,:)-ydot;
gradw=zeros(3);
gradtc=zeros(3,1);
gradts=zeros(3,1);

for k=1:3
    gradw(k,:)=sum([fit_error(:,1).*kerns(:,k) fit_error(:,2).*kerns(:,k) fit_error(:,3).*kerns(:,k)]); %After summing, should be a row vector

    ta=(tinds-tc(k))/ts+.5;
    tcdiff=(-60*ta+180ta.^2-120*ta.^3)/(ts(k).^2);
    gradtc(k)=sum(sum([fit_error(:,1).*tcdiff fit_error(:,2).*tcdiff fit_error(:,3).*tcdiff]));
    
    gradts(k)=sum(fit_error(:,1)*((-15/8)/(ts(k)^6)).*(2*t-2*tc(k)+ts(k)).*(-ts(k)+2*t-2*tc(k)).*(20*tc(k)^2-40*t*tc(k)-ts(k)^2+20*t.^2)*w(k,1));
end

cost=sum(sum(erroryfit.^2))+overlap;

grad=[gradw(:);gradtc;gradts];




%Finding points of low curvature and high velocity is equivalent to finding
%the corners of a puzzle?

%If you think you have a center and a direction worth starting from, in
%theory you just need a span. Any reason NOT to use change in dot product
%with time? Is there a way to test for symmetric rate of decline in
%direction instead BECAUSE you've presumed a stack depth?

%In theory, you can make this episodic via optimization 3-deep w,tc,ts and
%cost for making things more than 2-deep. Few free parameters.

%So sketch:
%1. Abuse our starting point with low curvature and use optimization in its
%tc+/-ts window to do gradient descent on w,tc,ts with the rule that the
%two other subs cannot overlap, which is just a logical*inf...ish.