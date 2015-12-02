clc
clear all

reachStruct=processPulse(1);

reachCats=[reachStruct.reachcat];

for RC=1 %unique(reachCats)
    subStruct=reachStruct(reachCats==RC);
    
    X=vertcat(subStruct.x);
    S=vecmag(vertcat(subStruct.v));
    
    x0=subStruct(1).x0;
    xf=subStruct(1).xf;
    
    if abs(xf(1)-x0(1))>.2
        rl=30;
    else
        rl=15;
    end
    X(:,1)=(X(:,1)-x0(1))/(xf(1)-x0(1))*(rl);
    X(:,2)=100*(X(:,2)-.5);
    
    xbins=linspace(0,rl,2*rl);
    ybins=linspace(-1,1,20);
    yhist=zeros(length(ybins)-1,length(xbins)-1);
    shist=yhist;
    for k=1:length(xbins)-1
        f=find((X(:,1)>xbins(k))&X(:,1)<(xbins(k+1)));
        for kk=1:length(ybins)-1
            ff=find((X(f,2)>ybins(kk))&X(f,2)<ybins(kk+1));
            yhist(kk,k)=length(ff);
            if length(ff)==0
                shist(kk,k)=0;
            else
                shist(kk,k)=mean(S(f(ff)));
            end
        end
        shist(:,k)=shist(:,k)/shist(length(ybins)/2-1,k);
    end
end

figure(1)
clf
surf(-log(yhist))
