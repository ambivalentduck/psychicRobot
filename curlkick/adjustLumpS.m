function x=adjustLumpS(x,S)

lL=length(S);

x(4*(1:lL)-3)=x(4*(1:lL)-3).*S./x(4*(1:lL)); %Lx
x(4*(1:lL)-2)=x(4*(1:lL)-2).*S./x(4*(1:lL)); %Ly
x(4*(1:lL))=S;