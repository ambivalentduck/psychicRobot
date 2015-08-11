function [shift,n,T]=fitShiftedGam(ts,varargin)

if nargin<2
    doPlotting=0;
else
    doPlotting=varargin{1};
end

global x

ts=ts(ts>0);
x=ts.^-2;

minx=min(x);
shiftEst=minx-.000001; %gamfit is weird about zeroes in the data, so use an epsilon
gpEst=gamfit(x-shiftEst);
gp=fminsearch(@shiftedGamObj,[shiftEst, gpEst]);
shift=gp(1);
n=gp(2);
T=gp(3);

if doPlotting
    xcdf=linspace(minx,max(x),25);
    gcdfEst=gamcdf(xcdf-shiftEst,gpEst(1),gpEst(2));
    gcdf=gamcdf(xcdf-shift,n,T);
    hold on
    ecdf(x,'bounds','on')
    plot(xcdf,gcdfEst,'m',xcdf,gcdf,'r')
    xlabel('ts^{-2}')
    ylabel('Cumulative Probability')
    title(['U=',num2str(shift),' n=',num2str(n),' T=',num2str(T)])
end

function cost=shiftedGamObj(params)

global x

cost=gamlike(params(2:3),x-params(1));
