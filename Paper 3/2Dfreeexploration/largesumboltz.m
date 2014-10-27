clc
clear all

Ws=[.35 2.085];

figure(1)
clf
N=10000;

sumvec=zeros(N,2);

for w=1:2
    W=Ws(w);
    invcdfboltz=@(x) -W*log(1-x)
    %sumvec(:,w)=cumsum(sign(rand(N,1)-.5).*invcdfboltz(.95*rand(N,1)));
    for t=2:N
        step=invcdfboltz(.95*rand);
        if sumvec(t-1,w)-step<0
            sumvec(t,w)=sumvec(t-1,w)+step;
        else
            sumvec(t,w)=sumvec(t-1,w)+sign(rand-.5)*step;
        end
    end
end

subplot(1,2,1)
hist((sumvec(:,1)),75)

subplot(1,2,2)
hist((sumvec(:,2)),75)


figure(2)
clf
t=0:.1:4;

plot(t,exp(-t.^20/.35))