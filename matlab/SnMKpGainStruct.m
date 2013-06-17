function out=SnMKpGainStruct(n)

if nargin<1
    out={1*[15 6; 6 16], 1.5*[15 6; 6 16], 3*[15 6; 6 16]};
else
    out=cell(length(n),1);
    for k=1:length(n)
        out{k}=n(k)*[15 6; 6 16];
    end
end