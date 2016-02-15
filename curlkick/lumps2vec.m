function x=lumps2vec(lumps)

x=zeros(length(lumps)*4,1);

for k=1:length(lumps)
    offset=(k-1)*4;
    x(offset+(1:2))=lumps(k).L;
    x(offset+3)=lumps(k).C;
    x(offset+4)=lumps(k).S;
end