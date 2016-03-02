function lumps=vec2lumps(x)

nlumps=(length(x)/4);

for k=nlumps:-1:1
    offset=(k-1)*4;
    lumps(k).L=x(offset+(1:2))';
    lumps(k).C=x(offset+3);
    lumps(k).S=x(offset+4);
end
