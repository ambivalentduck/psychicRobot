function [mue,n]=getMUE(bins,ref,y,method)

if nargin<4
    handle=@mean;
else
    switch lower(method)
        case 'max'
            handle=@max;
        otherwise
            handle=@mean;
    end
end

n=NaN*bins;

edges=(bins(2:end)+bins(1:end-1))/2;

x=y(:,1);

n(1)=handle(abs((y(x<=edges(1),2))-ref(1)));
for k=2:length(bins)-1
    n(k)=handle(abs(y((x<=bins(k))&(x>bins(k-1)),2)-ref(k)));
end
n(end)=handle(abs((y(x>edges(end),2)-ref(end))));

mue=handle(n(~isnan(n))); %max of maxes or mean of means, notice spacial reweighting if mean.