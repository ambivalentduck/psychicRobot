function correctcdf(trials)

% lambda=3;
M=10000;
% N=4;
% 
% c=.45*sum(sqrt(exprnd(lambda,M,N)),2);
% skewness(c)


tcats=[trials.targetcat];
dcats=[trials.disturbcat];

f=find((tcats==4)&~dcats);

segments=[.3 .2 .1 .05];
S=zeros(length(f),length(segments));
RW=S;

for k=1:length(f)
    tr=trials(f(k));
    [S(k,:),RW(k,:)]=getTsMetric(tr.x,tr.v,tr.a,tr.t,tr.x(1,:),tr.x(end,:),segments);
end

c=S(:,4);
c=sort(c);
c=c(1:end-1);
length(c)
skewness(c)
[F,x]=ecdf(c);

zed=linspace(0,max(x),M)';

    function P=computeP(params)
        lam=params(1);
        n=params(2);
        p=1/lam*exp(-zed.^2/lam)*2.*zed;
        p=real(ifft(fft(p).^n));
        P=cumsum(p);
        P=P/P(end);
        P=interp1(zed,P,x);
    end
minme=@(params) sum((F-computeP(params)).^2);

%minned=fmincon(minme,[.2; 2],-eye(2),[0;0])
minned=fmincon(minme,[.2; 2],[],[],[],[],[0;0],[],@(x) deal(0,0),optimoptions('fmincon','TolX',1e-14))
%minned=fminbnd(minme,[0 0],[1 100]) %,optimset('TolX',1e-12))

figure(1)
clf
hold on
plot(x,F,'b',x,computeP(minned),'r')
msePatton=mean((F-computeP(minned)).^2)
mu=mean(c);
sigma=std(c);
plot(x,normcdf(x,mu,sigma),'g')
mseGauss=mean((F-normcdf(x,mu,sigma)).^2)
mseGauss-msePatton

minned(1)

xlabel('Time, Seconds')
ylabel('Cumulative Probability')
legend('Data','Model','Gaussian')
end