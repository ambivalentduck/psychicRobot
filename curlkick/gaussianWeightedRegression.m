function [Y,SY]=gaussianWeightedRegression(x,y,X,K)

sy2=size(y,2);
Y=zeros(length(X),sy2);
lX=length(X);

for k=1:lX
    w=exp(-K*(x-X(k)).^2);
    sw=sum(w);
    for kk=1:sy2
        Y(k,kk)=sum(y(:,kk).*w)/sw;
        if nargout>1
            SY(k,kk)=sqrt(sum(w.*(y(:,kk)-Y(k,kk)).^2)/(sw*(lX-1)/lX));
        end
    end
end
