function C=lumps2rgbk(lumps,y)

%say that 1=r, 2=g 3=b, not those=black

sy1=size(y,1);

z=zeros(sy1,4);

for k=1:min(3,length(lumps)) %Counting to three is hazardous if you won't always succeed
    z(lumps(k).inds,4)=1;
    z(lumps(k).inds,k)=vecmag(lumps(k).y(:,3:4));
end

C=zeros(sy1,3);

f=find(z(:,4));
for kk=f'
    s=sum(z(kk,1:3));
    if s>0
        C(kk,:)=z(kk,1:3)/s;
    end
end
