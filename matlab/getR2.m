function out=getR2(X)

%PL=sum(sqrt(sum((X(2:end,:)-X(1:end-1,:)).^2,2)));
%out=PL/sqrt(sum((X(end,:)-X(1,:)).^2,2));

% Y=X(:,2);
% X=X(:,1);
% [b,trash1,trash2,trash2,stats]=regress(X+Y,[ones(length(X),1) X]);
% out=stats([1 3]);

out=[max(abs(X(:,2)-.5)) 1]*100;