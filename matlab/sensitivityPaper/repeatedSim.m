function simmed=repeatedSim(params,t,xvaf,handle,doplot)

global l1 l2 m1 m2 lc1 lc2 I1 I2 x0 kp0 kp1 kp kd kpgain kpkdratio reflexratio kpkdreflexratio

xvafnom=xvaf;

if doplot
    figure(4)
    clf
    hold on
end

sp1=size(params,1);
tic
for k=1:sp1
    v=params(k,:);
    l1=v(1);
    l2=v(2);
    lc1=l1*v(3);
    lc2=l2*v(4);
    mass=v(5);
    m1=v(6)*mass;
    m2=v(7)*mass;
    I1=m1*(v(8)*l1)^2;
    I2=m2*(v(9)*l2)^2;

    x0=[v(10) v(11)];
    xvaf(:,7:8)=[xvafnom(:,7)+v(12)+randn(length(t),1)*v(14) xvafnom(:,8)+v(13)+randn(length(t),1)*v(15)];

    kpgain=v(16);
    kp0=v(17)*[v(21) v(22); v(23) v(24)];
    kp1=v(18)*[v(25) v(26); v(27) v(28)];
    kp=v(19)*[v(29) v(30); v(31) v(32)];
    kd=v(20)*[v(33) v(34); v(35) v(36)];
    
    kpkdratio=v(37);
    reflexratio=v(38);
    kpkdreflexratio=v(39);

    simmed(N).y=extract(t,xvaf,handle);
    if doplot
        plot(simmed(N).y(:,1),simmed(N).y(:,2))
        axis equal
        drawnow
    end
    [k/sp1 toc/k ((sp1/k-1)*(toc))/60]
end