function params=paramsPopulator(varargin)

nom=[.33; %l1
    .34; %l2
    .436; %lc1/l1 factor, Winters 1990
    .682; %lc2/l2 factor, Winter 1990
    175*.4535; %mass, lb to kg conversion
    .028; %m1/mass conversion factor, Winters 1990
    .022; %m2/mass conversion factor, Winters 1990
    .322; %rog1 Winters 1990
    .468; %rog2 Winters 1990
    -0.0570; %Shoulder x in robot space
    .8800; %Shoulder y in robot space
    0; %Robot x bias
    0; %Robot y bias
    0; %Robot x noise
    0; %Robot y noise
    1; %Kp gain
    1; %Burdet const Kp gain
    1; %Burdet torque-varying Kp gain
    1; %Shad-muss Kp gain
    1; %Shad-muss Kd gain
    10.8; %Burdet kp0
    2.83;
    2.51;
    8.67;
    3.18; %Burdet kp1
    2.15;
    2.34;
    6.18;
    15; %Shad-Muss kp
    6;
    6;
    16;
    2.3; %Shad-Muss kd
    .09;
    .09;
    2.4
    1/12; %Kd/Kp for Burdet
    1/50; %Reflex Kp Ratio
    2]; % Kp/Kp ratio for reflexes
nom=nom';
lnom=length(nom);

sd=[.01; %l1, 1cm
    .01; %l2, 1cm
    .0695; %lc1/l1 factor, 1cm
    .0431; %lc2/l2 factor, 1cm
    3.1; %mass self-report error SD=3.1 kg (The Accuracy of Self-Reported Weights, stunkard and albaum Am J Clin Nutr 1981)
    .0029; %m1/mass conversion factor, Dempster 1955
    .00248; %m2/mass conversion factor, Dempster 1955
    .322/20; %rog1 5%
    .468/20; %rog2 5%
    .02; %Shoulder x in robot space
    .02; %Shoulder y in robot space
    .2310; %Robot x bias, empirical
    .1067; %Robot y bias, empirical
    .1653; %Robot x noise, empirical
    .2869; %Robot y noise, empirical
    .15; %Kp gain
    .15; %Burdet const Kp gain
    .15; %Burdet torque-varying Kp gain
    .15; %Shad-muss Kp gain
    .15; %Shad-muss Kd gain
    .15*10.8; %Burdet kp0
    .15*2.83;
    .15*2.51;
    .15*8.67;
    .15*3.18; %Burdet kp1
    .15*2.15;
    .15*2.34;
    .15*6.18;
    .15*15; %Shad-Muss kp
    .15*6;
    .15*6;
    .15*16;
    .15*2.3; %Shad-Muss kd
    .15*.09;
    .15*.09;
    .15*2.4
    .15*1/12; %Kd/Kp for Burdet
    .15*1/50; %Reflex Kp Ratio
    .15*2]; % Kp/Kp ratio for reflexes
sd=sd';

if nargin==3
    column=varargin{1};
    vals=varargin{2};
    lvals=length(vals);
    params=zeros(lvals,lnom);
    for k=1:lvals
        params(k,:)=nom;
    end
    switch lower(varargin{3})
        case 'multiply'
            params(:,column)=params(:,column).*vals;
        case 'add'
            params(:,column)=params(:,column)+vals;
        otherwise
            params(:,column)=vals;
    end
    return
end

n=varargin{1};

if length(n)==1
    span=[-3 -2 -1 -.5 0 .5 1 2 3];
    params=zeros(length(span),lnom);
    for k=1:length(span)
        params(k,:)=nom;
        params(k,n)=params(k,n)+span(k)*sd(n);
    end
    return
end

in=[1 1; %l1
    1 1; %l2
    1 1; %lc1/l1 factor
    1 1; %lc2/l2 factor
    1 1; %mass, lb to kg conversion
    1 1; %m1/mass conversion factor
    1 1; %m2/mass conversion factor
    1 1; %rog1
    1 1; %rog2
    1 1; %Shoulder x in robot space
    1 1; %Shoulder y in robot space
    1 1; %Robot x bias
    1 1; %Robot y bias
    1 1; %Robot x noise
    1 1; %Robot y noise
    0 0; %Kp gain
    0 1; %Burdet const Kp gain
    0 1; %Burdet torque-varying Kp gain
    1 0; %Shad-muss Kp gain
    1 0; %Shad-muss Kd gain
    0 0; %Burdet kp0
    0 0;
    0 0;
    0 0;
    0 0; %Burdet kp1
    0 0;
    0 0;
    0 0;
    0 0; %Shad-Muss kp
    0 0;
    0 0;
    0 0;
    0 0; %Shad-Muss kd
    0 0;
    0 0;
    0 0;
    0 1; %Kd/Kp for Burdet
    0 1; %Reflex Kp Ratio
    0 1]; % Kp/Kp ratio for reflexes

sumin=sum(in); % ==[17 20]

sn1=size(n,1);
sn2=size(n,2);

f=find(sn2==sumin);
if isempty(f)
    error('You''re misusing this function.')
end

params=zeros(sn1,lnom);

f2=find(in(:,f));

for k=1:length(lnom)
    if in(k,f) %Varying
        params(:,k)=norminv(n(:,f2==k),nom(k),sd(k));
    else
        params(:,k)=nom(k); %Set to nominal value
    end
end





