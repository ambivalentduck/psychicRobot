function [shift,n,T,rmse]=fitShiftedGam(x,varargin)

if isempty(x)
    error('Dude, x is emptier than a politician''s campaign promises.')
end

if nargin<2
    doPlotting=0;
else
    doPlotting=varargin{1};
end

minx=min(x);
shiftEst=minx-.000001; %gamfit is weird about zeroes in the data, so use an epsilon
gpEst=gamfit(x-shiftEst);
%gp=fmincon(@shiftedGamObj,[shiftEst gpEst]',diag([-1 -1 -1]),zeros(3,1));
gp=fminunc(@shiftedGamObj,[shiftEst gpEst]',optimset('TolX',1e-9,'MaxFunEvals',2000,'TolFun',1e-9));
shift=gp(1);
n=gp(2);
T=gp(3);

[f,xb]=ecdf(x);
gcdf=gamcdf(xb-shift,n,T);
rmse=mean(sqrt((f-gcdf).^2));

if doPlotting
    xcdf=linspace(minx,max(x),25);
    gcdfEst=gamcdf(xcdf-shiftEst,gpEst(1),gpEst(2));
    gcdf=gamcdf(xcdf-shift,n,T);
    hold on
    ecdf(x,'bounds','on')
    %plot(xcdf,gcdfEst,'m')
    plot(xcdf,gcdf,'r')
    ylabel('Cumulative Probability')
    title(['U=',num2str(shift),' n=',num2str(n),' T=',num2str(T)])
end

%Nesting the function allows sharing of x without using a global
    function cost=shiftedGamObj(params)
        cost=gamlike(params(2:3),x-params(1))+params(2)/2;
    end

end
