function [lump,y]=rulesFindLumps(t,y,ind)

%Rule 1:
% Submovements touch end-to-end in time. This guarantees a stack depth of two or less, but not necessarily pure peaks.

%Rule 2:
% Submovements are perfect, straight-line 5th order polynomials.

lT=length(t);
gT=t(2)-t(1); %Poor man's gradient, pray it's evenly sampled if you keep this.

dr=y(ind,3:4);
dr=dr/norm(dr);
dp1=zeros(lT,1);
dpunit1=dp1;
for kk=1:lT
    dp1(kk)=dot(y(kk,3:4),dr);
    dpunit1(kk)=dp1(kk)/norm(y(kk,3:4));
end
gu1=gradient(dpunit1);


minpeakheight=.01;
[vals,locs1]=findpeaks(gu1,'minpeakheight',minpeakheight);
[vals,locs2]=findpeaks(-gu1,'minpeakheight',minpeakheight);
lower1=locs1(find(locs1<ind,1,'last'));
upper1=locs2(find(locs2>ind,1,'first'));

if isempty(lower1)
    lower1=1;
end
if isempty(upper1)
    upper1=lT;
end

ispan=min(ind-lower1,upper1-ind);
inds=ind-ispan:ind+ispan;
tspan=t(inds(end))-t(inds(1));

% vmax=dist*1.875/tspan, see vc below.
span=y(ind,3:4)*tspan/1.875;
xi=y(ind,1:2)-span/2;
xf=y(ind,1:2)+span/2;

ta=t(inds)';
ta=ta-ta(1);
ta=ta/ta(end);
xc=(xi'*ones(size(ta))+(xf-xi)'*(10*ta.^3-15*ta.^4+6*ta.^5))';
vc=((xf-xi)'*(30*ta.^2-60*ta.^3+30*ta.^4)/tspan)';

lump.y=[xc vc];
lump.inds=inds;
lump.t=t(inds);

y(inds,:)=y(inds,:)-lump.y;

