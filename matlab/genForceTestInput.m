clc
clear all

KICKMAG=15;

%1 trial number
%2 x
%3 y
%4 early pulse mag and +/- direction
%5 late pulse mag and +/- direction
%6 pseudorandom white noise magnitude
%7 shape (-1=none,0=triangle,1=square, 2=circle, 3=inf desired)
%8 cursor shown (0/1 false/true)
%9 acquisitions needed



%10 reaches unperturbed
for k=1:10
    if mod(k,2)
        x=.15;
    else
        x=-.15;
    end
    out(k).dat=[x; .5; 0; 0; 0; -1; 1; 1];
end

%20 reaches kicked
for k=11:30
    if mod(k,2)
        x=.15;
    else
        x=-.15;
    end
    if rand>.5
        kick=KICKMAG;
    else
        kick=-KICKMAG;
    end
    if rand>.5
        early=kick;
        late=0;
    else
        early=0;
        late=kick;
    end
    out(k).dat=[x; .5; early; late; 0; -1; 1; 1];
end

%20 reaches white
for k=31:50
    if mod(k,2)
        x=.15;
    else
        x=-.15;
    end
    out(k).dat=[x; .5; 0; 0; 1.5; -1; 1; 1];
end

for k=0:3
    out(k+51).dat=[0; 0; 0; 0; 0; k; 1; 4];
end
for k=0:3
    out(k+55).dat=[0; 0; 0; 0; 0; k; 0; 4];
end
for k=0:3
    out(k+59).dat=[0; 0; 0; 0; 1; k; 0; 4];
end
for k=0:3
    out(k+63).dat=[0; 0; 0; 0; 1.5; k; 0; 4];
end
    

o=[out.dat]';
o=[(1:length(out))' o]
fid=fopen('../Data/input.dat','w');
fprintf(fid,'%5.0f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%1.0f\t%1.0f\t%1.0f\n',o');
fclose(fid);