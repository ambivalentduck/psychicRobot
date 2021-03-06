function params=paramsPopulator(varargin)

nom=[.353; %l1
    .363; %l2
    .436; %lc1/l1 factor, Winters 1990
    .682; %lc2/l2 factor, Winter 1990
    195*.4535; %mass, lb to kg conversion
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

if nargin==0
    params=nom;
    return
end

lnom=length(nom);

sd=[.017; %l1, Dempster 1955
    .011; %l2, Dempster 1955
    .15*.436; %lc1/l1 factor, 1cm
    .15*.682; %lc2/l2 factor, 1cm
    3.1; %mass self-report error SD=3.1 kg (The Accuracy of Self-Reported Weights, stunkard and albaum Am J Clin Nutr 1981)
    .0029; %m1/mass conversion factor, Dempster 1955
    .00248; %m2/mass conversion factor, Dempster 1955
    .15*.322; %rog1 5%
    .15*.468; %rog2 5%
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

names={'l1, m';
    'l2, m';
    'lc1/l1, unitless';
    'lc2/l2, unitless';
    'gross body mass, kg';
    'm1/mass, unitless';
    'm2/mass, unitless';
    'Seg1 RoG/l1, unitless';
    'Seg2 RoG/l2, unitless';
    'Shoulder x in robot space, m'
    'Shoulder y in robot space, m'
    'Force sensor bias x-axis, N';
    'Force sensor bias y-axis, N';
    'Force sensor Gaussian noise SD x-axis, N';
    'Force sensor Gaussian noise SD y-axis, N';
    'Kp gain, unitless';
    'Burdet const Kp gain, unitless';
    'Burdet torque-varying Kp gain, unitless';
    'Shad-Muss Kp gain, unitless';
    'Shad-Muss Kd gain, unitless';
    'Burdet kp0(1,1), N*M/rad';
    'Burdet kp0(1,2), N*M/rad';
    'Burdet kp0(2,1), N*M/rad';
    'Burdet kp0(2,2), N*M/rad';
    'Burdet kp1(1,1), N*M/rad';
    'Burdet kp1(1,2), N*M/rad';
    'Burdet kp1(2,1), N*M/rad';
    'Burdet kp1(2,2), N*M/rad';
    'Shad-Muss kp(1,1), N*M/rad';
    'Shad-Muss kp(1,2), N*M/rad';
    'Shad-Muss kp(2,1), N*M/rad';
    'Shad-Muss kp(2,2), N*M/rad';
    'Shad-Muss kd(1,1), N*M/rad';
    'Shad-Muss kd(1,2), N*M/rad';
    'Shad-Muss kd(2,1), N*M/rad';
    'Shad-Muss kd(2,2), N*M/rad';
    'Kd/Kp for Burdet, unitless';
    'Reflex Kp Ratio, unitless';
    'Kd/Kp ratio for reflexes, unitless'};

latexnames={'Upper Arm Length $(L_1)$';
    'Forearm Length $(L_2)$';
    'Upper Arm Center of Mass Ratio $(\frac{L_{m1}}{L_1})$';
    'Forearm Center of Mass Ratio $(\frac{L_{m2}}{L_2})$';
    'Gross Body Mass $(m_g)$';
    'Upper Arm Mass Ratio $(\frac{m_1}{m_g})$';
    'Forearm Mass Ratio $(\frac{m_2}{m_g})$';
    'Upper Arm Radius of Gyration Ratio $(\frac{ROG_1}{L_1})$';
    'Forearm Radius of Gyration Ratio $(\frac{ROG_2}{L_2})$';
    'Shoulder Parallel Coordinate'
    'Shoulder Perpendicular Coordinate'
    'Force Sensor Miscalibration $x$-axis';
    'Force Sensor Miscalibration $y$-axis';
    'Force Sensor Gaussian Noise SD $x$-axis';
    'Force Sensor Gaussian Noise SD $y$-axis';
    'Whole $K_P$ gain, unitless';
    'Torque-Invariant Impedance Misestimation Ratio';
    'Torque-Dependent Impedance Misestimation Ratio';
    'Shad-Muss Kp gain, unitless';
    'Shad-Muss Kd gain, unitless';
    '$K_{P0,11}, \frac{N\cdot M}{rad}$';
    '$K_{P0,12}, \frac{N\cdot M}{rad}$';
    '$K_{P0,21}, \frac{N\cdot M}{rad}$';
    '$K_{P0,22}, \frac{N\cdot M}{rad}$';
    '$K_{P1,11}, rad^{-1}$';
    '$K_{P1,12}, rad^{-1}$';
    '$K_{P1,21}, rad^{-1}$';
    '$K_{P1,22}, rad^{-1}$';
    'Shad-Muss kp(1,1), N*M/rad';
    'Shad-Muss kp(1,2), N*M/rad';
    'Shad-Muss kp(2,1), N*M/rad';
    'Shad-Muss kp(2,2), N*M/rad';
    'Shad-Muss kd(1,1), N*M/rad';
    'Shad-Muss kd(1,2), N*M/rad';
    'Shad-Muss kd(2,1), N*M/rad';
    'Shad-Muss kd(2,2), N*M/rad';
    'Damping-to-Stiffness Ratio $(k_d)$';
    'Reflex Impedance Scale Factor $(\frac{G}{K})$';
    'Reflex Damping-to-Stiffness Ratio $(g_d)$'};

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

in=[1 1 1 1; %l1
    1 1 1 1; %l2
    1 1 1 1; %lc1/l1 factor
    1 1 1 1; %lc2/l2 factor
    1 1 1 1; %mass, lb to kg conversion
    1 1 1 1; %m1/mass conversion factor
    1 1 1 1; %m2/mass conversion factor
    1 1 1 1; %rog1
    1 1 1 1; %rog2
    1 1 1 1; %Shoulder x in robot space
    1 1 1 1; %Shoulder y in robot space
    1 1 1 1; %Robot x bias
    1 1 1 1; %Robot y bias
    1 1 1 1; %Robot x noise
    1 1 1 1; %Robot y noise
    1 0 1 0; %Kp gain
    0 0 1 1; %Burdet const Kp gain
    0 0 1 1; %Burdet torque-varying Kp gain
    1 1 0 0; %Shad-muss Kp gain
    1 1 0 0; %Shad-muss Kd gain
    0 0 1 0; %Burdet kp0
    0 0 1 0;
    0 0 1 0;
    0 0 1 0;
    0 0 1 0; %Burdet kp1
    0 0 1 0;
    0 0 1 0;
    0 0 1 0;
    1 0 0 0; %Shad-Muss kp
    1 0 0 0;
    1 0 0 0;
    1 0 0 0;
    1 0 0 0; %Shad-Muss kd
    1 0 0 0;
    1 0 0 0;
    1 0 0 0;
    0 0 1 1; %Kd/Kp for Burdet
    0 0 1 1; %Reflex Kp Ratio
    0 0 1 1]; % Kp/Kp ratio for reflexes

sumin=sum(in); % ==[17 20]

if ischar(n)
    if strcmpi(n,'names')
        params=names;
        return
    end
    if strcmpi(n,'latex')
        params=latexnames;
        return
    end
    f=find(strcmpi(n,{'shadmuss','burdet'}));
    params=[nom' sd' in(:,2*f-1:2*f)];
    return
end

sn1=size(n,1);
sn2=size(n,2);

f=find(sn2==sumin);
if isempty(f)
    error('You''re misusing this function.')
end

params=zeros(sn1,lnom);

f2=find(in(:,f));

for k=1:lnom
    if in(k,f) %Varying
        params(:,k)=norminv(n(:,f2==k),nom(k),sd(k));
    else
        params(:,k)=nom(k); %Set to nominal value
    end
end





