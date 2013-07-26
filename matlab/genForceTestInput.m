clc
clear all

%x
%y
%early pulse mag and +/- direction
%late pulse mag and +/- direction
%pseudorandom white noise magnitude
%shape (-1=none,0=triangle,1=square, 2=circle, 3=inf desired)
%cursor shown (0/1 false/true)
%acquisitions needed

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
        kick=30;
    else
        kick=-30;
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
    out(k).dat=[x; .5; 0; 0; rand; -1; 1; 1];
end
    

o=[out.dat]';
o=[(1:length(out))' o];
fid=fopen('../Data/input.dat','w');
fprintf(fid,'%5.0f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%1.0f\t%1.0f\t%1.0f\n',o');
fclose(fid);