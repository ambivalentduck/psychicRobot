function R2=linfitR2(x,y)

inds=find(~(0==(1-y)));
lrf=log(1-y(inds));
x=x(inds);

mdl = fitlm(x,lrf);
R2=mdl.Rsquared.ordinary;
m=[x ones(size(x))]\lrf;
X=[min(x) max(x)];
plot(X,m(1)*X+m(2),'k')