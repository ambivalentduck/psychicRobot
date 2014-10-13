function [cost, grad]=sparseupMJPfit(P)

global xdot t Pknown

%% Step 0 P->w,ts

lP5=length(P)/5;

%This three defacto modes:
%1) lP5=2, two fit+one already known
%2) lP5=3, three fit
%3) lP5>3, trying to fit something like a whole movement at once.

r=reshape(P,5,lP5);

w=r(1:3,:);
tc=r(4,:);
ts=abs(r(5,:));

%% Step 1 cost....sum error.^2 over tc(1)+/-ts(1)

uppers=tc+ts/2;
lowers=tc-ts/2;

inds=find((t>lowers(1))&(t<uppers(1)));
tinds=t(inds);

ydot=zeros(length(inds),lP5);
kerns=zeros(length(inds),lP5);

for k=1:lP5
    ta=(tinds-tc(k))/ts(k)+.5;
    kerns(:,k)=(30*ta.^2-60*ta.^3+30*ta.^4)/ts(k);
    kerns(tinds<=lowers(k),k)=0;
    kerns(tinds>=uppers(k),k)=0;
    ydot=ydot+kerns(:,k)*(w(:,k)');
end

if exist('Pknown','var')
   %Contribute to ydot and calc overlap from this to 2 
else
    overlap=max(0, min(uppers(2),uppers(3))-max(lowers(2),lowers(3)));
    overlapdir=sign(tc(2)-tc(3));
end

fit_error=ydot-xdot(inds,:);
gradw=zeros(lP5,3);
gradtc=zeros(lP5,1);
gradts=zeros(lP5,1);

for k=1:lP5
    gradw(k,:)=sum([fit_error(:,1).*kerns(:,k) fit_error(:,2).*kerns(:,k) fit_error(:,3).*kerns(:,k)]); %After summing, should be a row vector

    ta=(tinds-tc(k))/ts(k)+.5;
    tcdiff=(-60*ta+180*ta.^2-120*ta.^3)/(ts(k)^2);
    gradtc(k)=sum(sum([w(1,k)*fit_error(:,1).*tcdiff, w(2,k)*fit_error(:,2).*tcdiff, w(3,k)*fit_error(:,3).*tcdiff]));
    
    td=(tinds-tc(k))/(ts(k)^2);
    tsdiff=(td.*(-60*ta+180*ta.^2-120*ta.^3)-kerns(:,k))/ts(k);
    gradts(k)=sum(sum([w(1,k)*fit_error(:,1).*tsdiff, w(2,k)*fit_error(:,2).*tsdiff, w(3,k)*fit_error(:,3).*tsdiff]));
%     if k~=1
%         gradts(k)=gradts(k)+overlap;
%     end
end

cost=sum(sum(fit_error.^2))/2+1*(overlap>0);

grad=zeros(5,lP5);
for k=1:lP5
    grad(:,k)=[gradw(k,:)';gradtc(k);gradts(k)];
end
grad=grad(:);