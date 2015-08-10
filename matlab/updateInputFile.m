clc
clear all

in=load('../Data/Data_pulse_old/input324.dat');

fid=fopen('updatedInput.dat','w');
fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\n',[in ones(size(in,1),1)]')
fclose(fid);
