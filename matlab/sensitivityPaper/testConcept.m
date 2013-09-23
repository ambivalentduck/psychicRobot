clc
clear all

N=10000;

proofofconcept=@(a,b,c) a;

p=sobolset(6,'Skip',1e3,'Leap',1e2); %double wide is necessary, rest are generic values to deal with idealty
p=scramble(p,'MatousekAffineOwen'); %Same. Cleans up some issues quickly and quietly

varied=p(1:N,:); %Generate a sobol-distributed [0-1] set that theoretically spans the space very very well.

varied=exp(-(varied-.5));

A=zeros(N,1);
B=A;
AB=A;

for k=1:N
    A(k)=proofofconcept(varied(k,1),varied(k,2),varied(k,3));
    B(k)=proofofconcept(varied(k,4),varied(k,5),varied(k,6));
    for V=1:3
        blah=varied(k,1:3);
        blah(V)=varied(k,V+3);
        AB(k,V)=proofofconcept(blah(1),blah(2),blah(3));
    end
end

for k=1:3
    E(k)=1/(2*N)*sum((A-AB(:,k)).^2);
    VA(k)=1/N*sum(B.*(AB(:,k)-A));
end

E
VA