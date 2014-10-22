function out=slfop(t,tc,ts)

lower=tc-ts/2;
upper=tc+ts/2;

ta=(t-lower)/ts;
out=[10*ta.^3-15*ta.^4+6*ta.^5,
    (30*ta.^2-60*ta.^3+30*ta.^4)/ts,
    (60*ta-180*ta.^2+120*ta.^3)/(ts^2)];

out(t<=lower,:)=0;
out(t>=upper,:)=1;


