function [Ut,Ux,Xgrid]=makeVPF(t,X,V,A)

%Make a centimeter grid for each arm so that discretization is easy: just
%subtract off the lower bound and round to the nearest centimeter

%Do this for each arm, maybe as a function

%know X, V, A as time series T-by-3

%Compute u(t)
udt=dot(V',A');
Ut=cumtrapz(t,udt);

%Blame u(t) on x(t) to infer u(x)

lower=min(X);
upper=max(X);
range=upper-lower;

binsize=.01*ones(1,3);

Xs=(0:binsize(1):range(1))';
Ys=(0:binsize(2):range(2))';
Zs=(0:binsize(3):range(3))';

Ux=cell(length(Xs),length(Ys),length(Zs));
Xgrid=Ux;

for x=1:length(Xs)
    for y=1:length(Ys)
        for z=1:length(Zs)
            Xgrid{x,y,z}=[Xs(x) Ys(y) Zs(z)]+lower;
        end
    end
end

for k=1:length(t)
    inds=floor((X(k,:)-lower)./binsize)+1;
    Ux{inds(1),inds(2),inds(3)}(end+1)=Ut(k);
end

end