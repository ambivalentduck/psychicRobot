function [C,owned]=lumps2rgbk(lumps)

%say that 1=r, 2=g 3=b, not those=black

% r g b c(gb) m(rb) y(rg)
cc=colorScheme();

scc1=size(cc,1);
ly=length(lumps(1).ownership);

z=zeros(ly,scc1+1);
for k=1:min(scc1,length(lumps)) %Counting to three is hazardous if you won't always succeed
    z(lumps(k).inds,end)=1;
    z(:,k)=lumps(k).ownership;
end

f=find(z(:,scc1+1));
C=zeros(ly,3);
for k=f'
    s=sum(z(k,1:scc1));
    if s==0
        continue
    end
    for kk=1:scc1
        C(k,:)=C(k,:)+cc(kk,:)*z(k,kk)/s;
    end
end

C(C>1)=1;
C(C<0)=0;
C(isnan(C))=0;

owned=zeros(ly,1);
owned(f)=1;