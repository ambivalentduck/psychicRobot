clc
clear all

N=2;
mode(1).suffix='';
mode(2).suffix='nograv';
mode(1).color='r';
mode(2).color='b';
mode(1).fitcolor=[1 .5 .5];
mode(2).fitcolor=[0 1 1];


figure(1)
clf
hold on

for M=1:length(mode)
    d=load(['3dfreeexp',num2str(N),mode(M).suffix,'.dat']);
    
    t=d(:,1);
    x=d(:,2:4);
    
    v=zeros(size(x));
    gT=gradient(t);
    for k=1:3
        v(:,k)=gradient(x(:,k))./gT;
    end
    
    bounded=(x(:,3)<.75)&(x(:,3)>.35);
    v=v(bounded,:);
    x=x(bounded,:);
    
    %[f,h]=ecdf(x(:,3));
    h=sort(x(:,3));
    f=(1:length(h))'/length(h);
    
    handles(M)=plot(h,log(1-f),[mode(M).color,'-'],'markersize',.01);
    
    fitinds=(f~=1)&(h<.6);
    h=h(fitinds);
    f=f(fitinds);
    X=[h ones(size(h))];
    Y=log(1-f);
    p=X\Y;
    yhat=X*p;
    EY=mean(Y);
    R2=1-sum((Y-yhat).^2)/sum((Y-EY).^2)
    plot(X(:,1),yhat,'color',mode(M).fitcolor,'linewidth',2)
    
end
legend(handles,{'No Support','Arm Support'})
xlabel('Height')
ylabel('ln Prob(height)')
