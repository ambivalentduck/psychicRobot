function [kern,progress]=MJ5P(t,tc,ts)

if size(t,2)>size(t,1)
    t=t';
end

if size(tc,2)>size(tc,1)
    tc=tc';
end

ts2=ts/2;

td=t-tc;
ta=(td+ts2)/ts;

logind=(td>=-ts2)&(td<=ts2);
ta=ta(logind);
kern=zeros(length(t),1);
kern(logind)=(30*ta.^2-60*ta.^3+30*ta.^4)/ts;

if nargout>1
    progress=zeros(length(t),1);
    progress(logind,:)=10*ta.^3-15*ta.^4+6*ta.^5;
    progress(td>=ts2,1)=1;
end




