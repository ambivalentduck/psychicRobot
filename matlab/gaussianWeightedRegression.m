function Y=gaussianWeightedRegression(x,y,X,K)

sy2=size(y,2);
Y=zeros(length(X),sy2);

for k=1:length(X)
    w=exp(-K*(x-X(k)).^2);
    sw=sum(w);
    for kk=1:sy2
        Y(k,kk)=sum(y(:,kk).*w)/sw;
    end
end
