function fitOverlapTheta

we=.3228;
ne=1.1196;
Ae=0.0899;
mass=.5;
N=200;

M=5;
costs=zeros(M,1)
inip=[.1447 .4396];
bestp=inip;
bestcost=doSim(inip)
for kkk=1:M
    params{kkk}=fmincon(@doSim,bestp',-eye(2),zeros(2,1));
    costs(kkk)=doSim(params{kkk});
    if costs(kkk)<bestcost
        bestcost=costs(kkk);
        bestp=params{kkk};
    end
end
[v,i]=min(costs);
inicost=doSim(inip)
costs(i)
params{i}

    function cost=doSim(P)
        A=P(1);
        xp=P(2);
        
        L2Tn2=exprnd(1,N,1);
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
        
        etotal=zeros(size(t));
        ptotal=zeros(size(t));
        
        for k=1:N
            inds=find((t>=t1(k))&(t<=t2(k)));
            tc=(t1(k)+t2(k))/2;
            ta=(t(inds)-tc)+.5;
            kern=dx(k)*(30*ta.^2-60*ta.^3+30*ta.^4);
            acc=dx(k)*(60*ta-180*ta.^2+120*ta.^3);
            
            % Energy superposition
            ptotal(inds)=ptotal(inds)+mass*kern.*acc;
            etotal(inds)=etotal(inds)+.5*mass*kern.^2;
        end
        
        w=expfit(abs(ptotal));
        p=gamfit(etotal);
        n=p(1);
        A=p(2);
        cost=sum(abs([w n A]-[we ne Ae])./[we ne Ae]);
    end
end