function [kern,fop]=slmj5op(t,tc,ts)

if size(t,2)>size(t,1)
    t=t';
end

if size(tc,2)>size(tc,1)
    tc=tc';
end

ts2=ts/2;

td=t-tc;
ta=(td+ts2)/ts;

kern=(30*ta.^2-60*ta.^3+30*ta.^4)/ts;
kern(td<=-ts2)=0;
kern(td>=ts2)=0;

if nargout>1
    fop=[10*ta.^3-15*ta.^4+6*ta.^5,kern,(60*ta-180*ta.^2+120*ta.^3)/(ts^2)];
    fop(td<=-ts2,:)=0;
    fop(td>=ts2,1)=1;
    fop(td>=ts2,2:3)=0;
end


