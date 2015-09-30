function fitOverlapTheta2D

we=.3228;
ne=1.1196;
Ae=0.0899;
N=200;
L2Tn2=exprnd(1,N,1);

M=1;
costs=zeros(M,1)
inip=[.5706];
inicost=doSim(inip)
for kkk=1:M
    params{kkk}=fmincon(@doSim,inip',-eye(1),zeros(1,1));
    costs(kkk)=doSim(params{kkk});
end
[v,i]=min(costs);
inicost=doSim(inip)
costs(i)
params{i}

    function cost=doSim(P)
        A=.1396;
        xp=.447;
        mass=.5;
        dTheta=P(1);
        
        dx=sqrt(L2Tn2*A);
        
        t1=zeros(N,1);
        t2=t1;
        for k=1:N
            t2(k)=t1(k)+1;
            
            if k~=N
                t1(k+1)=t1(k)+xp;
            end
        end
        t=0:.005:t2(end);
        
        vtotal=zeros(length(t),2);
        atotal=zeros(length(t),2);
        
        lastTheta=0;
        for k=1:N
            inds=find((t>=t1(k))&(t<=t2(k)));
            tc=(t1(k)+t2(k))/2;
            ta=(t(inds)-tc)+.5;
            theta=dTheta*randn+lastTheta;
            lastTheta=theta;
            ux=[cos(theta); sin(theta)];
            kern=ux*dx(k)*(30*ta.^2-60*ta.^3+30*ta.^4);
            acc=ux*dx(k)*(60*ta-180*ta.^2+120*ta.^3);
            
            % Cartesian superposition
            vtotal(inds,:)=vtotal(inds,:)+kern';
            atotal(inds,:)=atotal(inds,:)+acc';
        end
        cartPm=mass*dot(atotal',vtotal')';
        cartEm=.5*mass*dot(vtotal',vtotal')';
        cartVm=vtotal;
        
        
        w=expfit(abs(cartPm));
        p=gamfit(cartEm);
        n=p(1);
        A=p(2);
        cost=sum(abs([w n A]-[we ne Ae])./[we ne Ae])
    end
end